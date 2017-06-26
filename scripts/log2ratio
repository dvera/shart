#!/usr/bin/env bash

usage() {
  echo "Usage: $(basename $0) [-t threads ] -e sample1_early.bg,[sample2_early.bg ... sampleN_early.bg ] -l sample1_late.bg,[sample2_late.bg ... sampleN_late.bg" 1>&2
  echo "" 1>&2
  exit 1
}

while getopts ":t:" opt; do
  case $opt in
  t)
    NTHREADS=$OPTARG
   ;;
  e)
    EBGS=$OPTARG
   ;;
  l)
    LBGS=$OPTARG
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

if [[ -z $EBGS ]] ; then
  echo 'no early files given'
  exit 1
fi
if [[ -z $LBGS ]] ; then
  echo 'no late files given'
  exit 1
fi

#INPUT=$@

l2r(){
  EBG=$1
  LBG=$2
  # check if file exists
  if [[ ! -f $EBG ]]; then echo "file \"$EBG\" not found"; exit 1; fi
  if [[ ! -f $LBG ]]; then echo "file \"$LBG\" not found"; exit 1; fi
  
  # check input is a fastq
  if [[ $EBG != *.bg ]]; then echo "input is not a bedGraph file"; exit 1; fi
  if [[ $LBG != *.bg ]]; then echo "input is not a bedGraph file"; exit 1; fi
  
  OUTPUT="$(basename $EBG | sed 's/\.bg$/_l2r\.bg/g')"
  LOGFILE=${OUTPUT}.log
  
  paste $EBG $LBG | awk '{print $1,$2,$3,log(($4+1)/($8+1))/log(2)} }' OFS='\t' > $OUTPUT
  
  echo $OUTPUT
}

export -f l2r

EBGS=$(echo $EBGS | tr ',' ' ')
LBGS=$(echo $LBGS | tr ',' ' ')

parallel --will-cite -k -j $NTHREADS "l2r {}" ::: $INPUT


# BGLINES=$(wc -l $EARLYBG | cut -d" " -f 1)
# RTLINES=$(wc -l ${PREFIX}.bg | cut -d" " -f 1)
# PCTZEROES=$(echo "${RTLINES}/${BGLINES}" | bc -l)
# echo "${PCTZEROES} of the ${BGLINES} windows have 0 reads in either fraction" > ${PREFIX}.bg.log
