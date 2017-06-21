#!/usr/bin/env bash

usage() {
  echo "Usage: $(basename $0) [-t threads ] fastq1 [fastq2 ... fastqN ]" 1>&2
  echo "" 1>&2
  exit 1
}

while getopts ":t:" opt; do
  case $opt in
  t)
    NTHREADS=$OPTARG
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

FASTQFILES=$@

# check if all files exist
for f in $FASTQFILES; do
  if [[ ! -f $f ]]; then
    echo "fastq file \"$f\" not found"
    exit 1
  fi
done


for f in $FASTQFILES; do
 OUTFILE="$(basename $f | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')_clip.fastq"
 LOGFILE=${OUTFILE}.log
 echo "cutadapt -a AGATCGGAAGAGCACACGTCTG -q 0 -O 1 -m 0 -o $OUTFILE $f > $LOGFILE && fastqc -q $OUTFILE"
done | parallel --citation -j $NTHREADS


 
 
