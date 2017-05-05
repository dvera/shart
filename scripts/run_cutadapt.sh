#!/usr/bin/env bash

INPUT=$1
PREFIX=$2

cutadapt -a AGATCGGAAGAGCACACGTCTG -q 0 -O 1 -m 0 -o ${PREFIX}_clip.fastq $INPUT > ${PREFIX}_clip.fastq.log fastqc $INPUT
fastqc ${PREFIX}_clip.fastq
