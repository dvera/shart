#!/usr/bin/env bash

EARLYBG=$1
LATEBG=$2
PREFIX=$3

paste $EARLYBG $LATEBG | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)} }' OFS='\t' > ${PREFIX}.bg
BGLINES=$(wc -l $EARLYBG | cut -d" " -f 1)
RTLINES=$(wc -l ${PREFIX}.bg | cut -d" " -f 1)
PCTZEROES=$(echo "${RTLINES}/${BGLINES}" | bc -l)
echo "${PCTZEROES} of the ${BGLINES} windows have 0 reads in either fraction" > ${PREFIX}.bg.log
