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
  echo 'no bam files specified'
  exit 1
fi

# DEBUG

FASTQFILES=$@

# check if all files exist
for f in $FASTQFILES; do
  if [[ ! -f $f ]] ; then
    echo "bam file \"$f\" not found"
    exit 1
  fi
done

# check if all files exist
for f in $FASTQFILES; do
  if [[ $f != *.bam ]]; then
    echo "file \"$f\" is not a bam file"
    exit 1
  fi
done

echo "removing duplicates"
for f in $FASTQFILES; do
  OUTPUT=$(basename $f | sed 's/\.bam$/_rmdup\.bam/g')
  echo "samtools rmdup -s $f $OUTPUT && samtools stats ${OUTPUT} > ${OUTPUT}.samstats"
done | parallel -j $NTHREADS

for f in $FASTQFILES; do
  OUTPUT=$(basename $f | sed 's/\.bam$/_rmdup\.bam/g')
  #echo "SAMLINES=$(samtools view -c $f);RSAMLINES=$(samtools view -c $OUTPUT);PCTNONZEROES=$(echo "${RSAMLINES}/${SAMLINES}" | bc -l)
  echo $(echo "$(samtools view -c $OUTPUT)/$(samtools view -c $f)" | bc -l) > ${OUTPUT}.pct
done | parallel -j $NTHREADS

  



 
 
