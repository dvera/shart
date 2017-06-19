#!/usr/bin/env bash

INPUT=$1
OUTPUT=$(basename $INPUT | sed 's/\.bam$/_rmdup\.bam/g')

samtools rmdup -s $INPUT $OUTPUT

samtools stats $OUTPUT > ${OUTPUT}.samstats

SAMLINES=$(samtools view -c $INPUT)

RSAMLINES=$(samtools view -c $OUTPUT)

PCTNONZEROES=$(echo "${RSAMLINES}/${SAMLINES}" | bc -l)

echo "${PCTNONZEROES} of the ${SAMLINES} reads were kept after removing duplicates" | tee ${OUTPUT}.log
