#!/bin/bash

usage() {
  echo "Usage: $(basename $0) [-t threads] -i bwaIndex fastq1 [fastq2 ... fastqN ]" 1>&2
  echo "" 1>&2
  echo "  bwaIndex can be a path to a bwa index prefix or a tarball of an bwa index" 1>&2
  echo "" 1>&2
  exit 1
}

while getopts ":i:t:" opt; do
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

# DEBUG
echo "index is $INDEXFILE"

FASTQFILES=$@

if [[ -f ${INDEXFILE}.bwt ]]; then
  INDEXPREFIX=$INDEXFILE
elif [[ -f $INDEXFILE ]] && [[ $INDEXFILE =~ "\.tar\.gz$" ]]; then
  mkdir -p bwaIndex
  INDEXPATH=$(readlink -f $INDEXFILE)
  tar -C bwaIndex -x -f $INDEXPATH
  INDEXPREFIX=$(readlink -f bwaIndex/*.bwt | sed 's/\.bwt$//g')
else
  echo "index not found"
  exit 2
fi

# run bwa
for f in $FASTQFILES; do
  # unzip file if compressed
  echo "processing $f"
  if [[ ! -f $f ]]; then
    echo "fastq file not found"
  exit 1
  elif [[ $f =~ \.gz$ ]]; then
    fo=${fastq/%\.gz/}
    gunzip -c $f > $fo
    f=$fo
  fi

  # create output prefix
  OUTPUT="$(basename $f | sed 's/\.fq$//g' | sed 's/\.fastq$//g').bam"
  
  # check to see input name has same as output name
  #  if [[ $OUTPUT == $f ]]; then
  #   echo "fastq file has inappropriate name, must be a fastq"
  #   exit 1
  #  fi
  
  # align fastq file and run samstats
  echo "running the following command:"
  echo "bwa mem -t $NTHREADS $INDEXPREFIX $f | samtools view -Shb - > $OUTPUT"
  bwa mem -t $NTHREADS $INDEXPREFIX $f | samtools view -Shb - > $OUTPUT
  samtools stats $OUTPUT > ${OUTPUT}.samstats
done
 
 
