# Docker-4dn-hic

This repo contains the source files for a docker image stored in duplexa/4dn-hic:v8\. (we will change the docker hub account soon)

## Table of contents

- [Cloning the repo](#cloning-the-repo)
- [Tool specifications](#tool-specifications)
- [Building docker image](#building-docker-image)
- [Sample data](#sample-data)
- [Tool wrappers](#tool-wrappers)

  - [run-list.sh](#run-listsh)
  - [run-bwa-mem-single.sh](#run-bwa-memsh)
  - [run-sort-bam.sh](#run-sort-bamsh)

## Cloning the repo

```
git clone https://github.com/4dn-dcic/docker-4dn-hic
cd docker-4dn-hic
```

## Sample data

Sample data files that can be used for testing the tools are included in the `sample_data` folder. These data are not included in the docker image.

## Tool wrappers

Tool wrappers are under the `scripts` directory and follow naming conventions `run-xx.sh`. These wrappers are copied to the docker image at built time and may be used as a single step in a workflow.

```
# default
docker run vera/docker-4dn-repliseq

# specific run command
docker run vera/docker-4dn-repliseq <run-xx.sh> <arg1> <arg2> ...

# may need -v option to mount data file/folder if they are used as arguments.
docker run -v /data1/:/d1/:rw -v /data2/:/d2/:rw vera/docker-4dn-repliseq <run-xx.sh> /d1/file1 /d2/file2 ...
```

### run-list.sh

Default command for this docker image. It lists the run commands available.

### run-bwa-mem-single.sh

Alignment module for single-end repli-seq data, based on bwa-mem.

- Input : a fastq file for an early- or late-enriched library
- Output : a bam file (Lossless, not sorted by coordinate)

#### Usage

Run the following in the container.

```
run-bwa-mem.sh <fastq1> <fastq2> <bwaIndex> <output_prefix> <nThreads>
# fastq1, fastq2 : input fastq files, either gzipped or not
# bwaIndex : tarball for bwa index, .tgz.
# output_prefix : prefix of the output bam file.
# nThreads : number of threads
```

# initial QC

run_fastqc.sh INPUT

# clip and reqc

run_cutadapt.sh INPUT PREFIX

# align, sort, get stats

run_bwa-mem-single.sh INPUT INDEX PREFIX NTHREADS

# remove duplicates and get stats

run_samtools-rmdup.sh INPUT PREFIX

# make windows

run_bedtools-makewindows.sh CHROMSIZES WINDOWSIZE PREFIX

# calcaulte coverage

run_bedtools-intersect.sh INPUT WINDOWS PREFIX

# calculate log2 ratios

run_log2ratio.sh EARLYBG LATEBG PREFIX

# create a reference distribution from multiple files

run_create-reference-distribution.sh PREFIX FILE1 FILE2 FILE3

# quantile normalize

run_quantile-normalization.sh INPUT REFERENCE PREFIX
