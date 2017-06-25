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

while getopts ":i:t:g:p" opt; do
  case $opt in
  e)
   EARLYBGS=$OPTARG
   ;;
  l)
   LATEBGS=$OPTARG
   ;;
  p)
   PAIRED=1
   ;;
  t)
   NTHREADS=$OPTARG
   ;;
  w)
   WINDOWSIZE=$OPTARG
   ;;
  c)
   CHROMSIZES=$OPTARG
   ;;
  i)
   INDEXFILE=$OPTARG
   ;;
  m)
   MEMPERTHREAD=$OPTARG
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

shift $((OPTIND-1))

if [[ -z $WINDOWSIZE ]]; then
  WINDOWSIZE=5000
fi

if [[ -z $MEMPERTHREAD ]]; then
  MEMPERTHREAD=5G
fi

if [[ -z $INDEXFILE ]]; then
  echo "must define an fasta or bwa index with -i"
fi

if [[ -z $CHROMSIZES ]]; then
  if [[ $INDEXFILE == *.fa* ]]; then
    awk '{
      if($0~">"){
        if(length(chrom)!=0){
          print chrom,clen
        }
        chrom=$1;
        clen=0
       } else{
        clen+=length($0)
       }
      } END{
        print chrom,clen
    }' OFS='\t' $INDEXFILE | tr -d '>' | sort -m -k1,1 > genome.chrom.sizes
    CHROMSIZES="genome.chrom.sizes"
  else
    echo "must define fasta as index if no chrom sizes file defined"
    exit 1
  fi
else
  awk 'if($2 ~ /^[0-9]+$/ || NF != 2) { exit 1 }'
  if [[ $? -gt 0 ]];
    echo "chromsizes looks incorrect"
    exit 1
   fi
fi

if [[ $CHROMSIZES == *.fa* ]]; then
    awk '{
      if($0~">"){
        if(length(chrom)!=0){
          print chrom,clen
        }
        chrom=$1;
        clen=0
       } else{
        clen+=length($0)
       }
      } END{
        print chrom,clen
    }' OFS='\t' $CHROMSIZES | tr -d '>' | sort -m -k1,1 > genome.chrom.sizes
    CHROMSIZES="genome.chrom.sizes"
fi

if [ -z $NTHREADS ]; then
  NTHREADS=1
fi

if [[ $# -eq 0 ]] ; then
  echo 'no fastq files specified'
  exit 1
fi

#check that dependencies are in PATH and correct versions
#check if chromsizes or fasta specified with -g
#check input files

if [[ -z $PAIRED ]]; then
  #todo: make sure length of EARLYBGS and LATEBGS is even
  #todo: check to see if read names in r1 and r2 are identical
  ER1=$(echo $EARLYBGS | tr ',' '\n' | paste - - | cut -f 1)
  ER2=$(echo $EARLYBGS | tr ',' '\n' | paste - - | cut -f 2)
  LR1=$(echo $LATEBGS | tr ',' '\n' | paste - - | cut -f 1)
  LR2=$(echo $LATEBGS | tr ',' '\n' | paste - - | cut -f 2)
else
  ER1=$(echo $EARLYBGS | tr ',' '\n')
  LR1=$(echo $LATEBGS | tr ',' '\n')
fi




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
for f in $ER1 $LR1 $ER2 $LR2; do
 BASE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')"
 OUTPUT=${BASE}_clip.fastq
 LOGFILE=${OUTPUT}.log
 echo "cutadapt -a AGATCGGAAGAGCACACGTCTG -q 0 -O 1 -m 0 -o $OUTPUT $f > $LOGFILE && fastqc -q $OUTPUT"
done | parallel --will-cite -j $NTHREADS

################################
### ALIGN READS WITH BWA #######
################################

# run bwa
for  in $FASTQFILES; do
  # create output prefix
  OUTPUT="$(basename $f | sed 's/\.fq$//g' | sed 's/\.fastq$//g').bam"
  # align fastq file and run samstats
  echo "bwa mem -v 2 -t $NTHREADS $INDEXPREFIX $f | samtools view -Shb - 2> ${OUTPUT}.log > $OUTPUT"
  bwa mem -v 2 -t $NTHREADS $INDEXPREFIX $f | samtools view -Shb - 2> ${OUTPUT}.log > $OUTPUT
done < 

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
  echo "samtools view -bq 20 $INPUT | samtools sort -m $MEMPERTHREAD -T $OUTPUT - > ${OUTPUT}"
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
  echo "bedtools bamtobed -i $INPUT | cut -f 1,2,3,4,5,6 | sort -S $MEMPERTHREAD -T . -k1,1 -k2,2n > $OUTPUT"
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

