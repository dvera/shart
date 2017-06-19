#!/usr/bin/env bash

INPUT=$1
RTBG=${@:3}
REFERENCE=$2

for INPUT in $RTBG; do
  OUTPUT=$(basename $INPUT | sed 's/\.bg/_qnorm\.bg/g')
  BGLINES=$(wc -l $INPUT | cut -d" " -f 1)
  paste <(sort -T . -k4,4g $INPUT) <(shuf -n $BGLINES $REFERENCE | sort -k4,4g) | \
    cut -f 1,2,3,8 | sort -k1,1 -k2,2n > $OUTPUT
done
