#!/usr/bin/env bash

INPUT=$1
PREFIX=$3

samtools rmdup -s $INPUT ${PREFIX}_rmdup.bam
samtools stats $INPUT > ${INPUT}.samstats
SAMLINES=$(wc -l $INPUT | cut -d" " -f 1)
RSAMLINES=$(wc -l ${PREFIX}_rmdup.bam | cut -d" " -f 1)
PCTNONZEROES=$(echo "${RSAMLINES}/${SAMLINES}" | bc -l)
echo "${PCTNONZEROES} of the ${SAMLINES} reads were kept after removing duplicates" > ${PREFIX}_rmdup.bam.log
