#!/usr/bin/env bash

INPUT=$1
PREFIX=$2

samtools rmdup -s $INPUT ${PREFIX}_rmdup.bam

samtools stats ${PREFIX}_rmdup.bam > ${PREFIX}_rmdup.bam.samstats

SAMLINES=$(samtools view -c $INPUT)

RSAMLINES=$(samtools view -c ${PREFIX}_rmdup.bam)

PCTNONZEROES=$(echo "${RSAMLINES}/${SAMLINES}" | bc -l)

echo "${PCTNONZEROES} of the ${SAMLINES} reads were kept after removing duplicates" | tee ${PREFIX}_rmdup.bam.log
