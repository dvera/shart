#!/bin/bash

usage() {
  echo "Usage: $(basename $0) [-t threads] -i bwaIndex fastq1 [fastq2 ... fastqN ]" 1>&2
  echo "" 1>&2
  echo "  bwaIndex can be a path to a bwa index prefix or a tarball of an bwa index" 1>&2
  echo "" 1>&2
  exit 1
}

################################
### PARSE COMMAND LINE ARGS ####
################################

while getopts ":i:t:g:" opt; do
  case $opt in
  t)
    NTHREADS=$OPTARG
   ;;
  i)
   INDEXFILE=$OPTARG
   ;;
  \?)
   echo "Invalid option: -$OPTARG" >&2
   usage
   ;;
  [?])
   usage
   ;;
  :)
   echo "Option -$OPTARG requires an argument." >&2
   echo "" >&2
   usage
   ;;
  esac
done

if [ -z $NTHREADS ]; then
  NTHREADS=1
fi

shift $((OPTIND-1))

if [[ $# -eq 0 ]] ; then
  echo 'no fastq files specified'
  exit 1
fi

#check if chromsizes or fasta specified with -g
#check input files

echo "index is $INDEXFILE"

FASTQFILES=$@

################################
### EXAMINE BWA INDEX ##########
################################

if [[ -f ${INDEXFILE}.bwt ]]; then
  INDEXPREFIX=$INDEXFILE
elif [[ -f $INDEXFILE ]] && [[ $INDEXFILE == *.tar.gz ]]; then
  mkdir -p bwaIndex
  INDEXPATH=$(readlink -f $INDEXFILE)
  tar -C bwaIndex -x -f $INDEXPATH
  INDEXPREFIX=$(readlink -f bwaIndex/*.bwt | sed 's/\.bwt$//g')
elif  [[ -f $INDEXFILE ]] && [[ $INDEXFILE == *.fa* ]]; then
  echo "fasta defined for index, making bwa index from fasta"
  mkdir -p bwaIndex
  INDEXFILE=$(readlink -f $INDEXFILE)
  (cd bwaIndex && bwa index -p genome $INDEXFILE)
  INDEXFILE=$(readlink -f bwaIndex)/genome
else
  echo "index not found"
  exit 2
fi

################################
### CLIP ADAPTERS FROM READS ###
################################

echo "clipping adapters from reads"
for f in $FASTQFILES; do
 BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
 OUTPUT=${BASE}_clip.fastq
 LOGFILE=${OUTPUT}.log
 echo "cutadapt -a AGATCGGAAGAGCACACGTCTG -q 0 -O 1 -m 0 -o $OUTPUT $f > $LOGFILE && fastqc -q $OUTPUT"
done | parallel --will-cite -j $NTHREADS

################################
### ALIGN READS WITH BWA #######
################################

# run bwa
for f in $FASTQFILES; do
  # create output prefix
  OUTPUT="$(basename $f | sed 's/\.fq$//g' | sed 's/\.fastq$//g').bam"
  # align fastq file and run samstats
  echo "bwa mem -v 2 -t $NTHREADS $INDEXPREFIX $f | samtools view -Shb - 2> ${OUTPUT}.log > $OUTPUT"
  bwa mem -v 2 -t $NTHREADS $INDEXPREFIX $f | samtools view -Shb - 2> ${OUTPUT}.log > $OUTPUT
done

echo "calculating alignment statistics"
for f in $FASTQFILES; do
  BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
  INPUT=${BASE}.bam
  OUTPUT=${BASE}.bam.samstats
  echo "samtools stats $INPUT > $OUTPUT"
done | parallel --will-cite -j $NTHREADS

################################
### SORT/FILTER ALIGNMENTS #####
################################

echo "sorting and filtering alignments"
for f in $FASTQFILES; do
  BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
  INPUT=${BASE}.bam
  OUTPUT=${BASE}_q20_sort.bam
  echo "samtools view -bq 20 $INPUT | samtools sort -T $OUTPUT - > ${OUTPUT}"
done | parallel --will-cite -j $NTHREADS

echo "calculating alignment statistics"
for f in $FASTQFILES; do
  BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
  INPUT=${BASE}_q20_sort.bam
  OUTPUT=${BASE}_q20_sort.bam.samstats
  echo "samtools stats $INPUT > $OUTPUT"
done | parallel --will-cite -j $NTHREADS

################################
### REMOVE PCR DUPLICATES ######
################################

echo "removing duplicates"
for f in $FASTQFILES; do
  BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
  INPUT=${BASE}_q20_sort.bam
  OUTPUT=${BASE}_q20_sort_rmdup.bam
  LOGFILE=${BASE}_q20_sort_rmdup.bam.pct
  echo "samtools rmdup -s $INPUT $OUTPUT 2> $LOGFILE"
done | parallel --will-cite -j $NTHREADS

echo "calculating alignment statistics"
for f in $FASTQFILES; do
  BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
  INPUT=${BASE}_q20_sort_rmdup.bam
  OUTPUT=${BASE}_q20_sort_rmdup.bam.samstats
  echo "samtools stats $INPUT > $OUTPUT"
done | parallel --will-cite -j $NTHREADS

################################
### CONVERT BAM TO BED #########
################################

echo "converting bam to bed"
for f in $FASTQFILES; do
  BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
  INPUT=${BASE}_q20_sort_rmdup.bam
  OUTPUT=${BASE}_q20_sort_rmdup.bed
  echo "bedtools bamtobed -i $INPUT | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n > $OUTPUT"
done | parallel --will-cite -j $NTHREADS


echo "calculating read densities"
for f in $FASTQFILES; do
  BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
  INPUT=${BASE}_q20_sort_rmdup.bed
  OUTPUT=${BASE}_q20_sort_rmdup.bg
  echo "bedtools intersect -sorted -c -b $INPUT -a $WINDOWBED | awk -v SCALE=$SCALE '{print $1,$2 ,$3 ,$4*SCALE }' OFS='\t'"
done | parallel --will-cite -j $NTHREADS


  # unzip file if compressed
#   echo "processing $f"
#   if [[ ! -f $f ]]; then
#     echo "fastq file not found"
#   exit 1
#   elif [[ $f == *.gz ]]; then
#     echo "extracting fastq files"
#     fo=${f%.*}
#     gunzip -c $f > $fo
#     f=$fo
#   fi

