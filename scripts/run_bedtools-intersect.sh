#!/usr/bin/env bash

INPUT=$1
WINDOWBED=$2
PREFIX=$3

bedtools bamtobed -i $INPUT | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n > ${PREFIX}.bed
SCALE=$(echo "1000000/$(wc -l $INPUT | cut -d" " -f1)" | bc -l)
bedtools intersect -sorted -c -b $INPUT -a $WINDOWBED | awk -v SCALE=$SCALE '{print $1,$2 ,$3 ,$4*SCALE }' OFS='\t' > ${PREFIX}.bg
