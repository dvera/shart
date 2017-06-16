#!/usr/bin/env bash

INPUT=$1
INDEX=$2
PREFIX=$3
NTHREADS=$4

run-bwa-mem-single.sh $INPUT $INDEX $PREFIX $NTHREADS run-sort-bam.sh $INPUT $PREFIX samtools stats $INPUT > ${INPUT}.samstats
