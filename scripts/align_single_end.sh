#!/bin/bash

fastq=$1
index_file=$2
prefix=$3
nThreads=$4

# unzip index
tar -xzf $index_file
index=`ls -1 *.bwt | head -1 | sed 's/\.bwt//g'`

# unzip fastq files
if [[ $fastq =~ \.gz$ ]]
then
  gunzip $fastq
  fastq=${fastq/%\.gz/}
fi

# run bwa
bwa mem -t $nThreads $index $fastq | samtools view -Shb - > $prefix.bam

samtools stats ${prefix}.bam > ${prefix}.bam.samstats
