#!/usr/bin/env bash

INPUT=$1
REFERENCE=$2
PREFIX=$3

BGLINES=$(wc -l $INPUT | cut -d" " -f 1)
paste <(sort -T . -k4,4g $INPUT) <(shuf -n $BGLINES $REFERENCE | sort -k4,4g) | cut -f 1,2,3,8 | sort -T . -k1,1 -k2,2n > ${PREFIX}_qnorm.bg
