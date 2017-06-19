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
cd docker-4dn-repliseq
```

## Sample data

Sample data files that can be used for testing the tools are included in the `sample_data` folder. These data are not included in the docker image.

## Tool wrappers

Tool wrappers are under the `scripts` directory and follow naming conventions `run_xx.sh`. These wrappers are copied to the docker image at build time and may be used as a single step in a workflow.

# Usage

### example
```bash
# execute a step on data in the current directory
docker run -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq <run-xx.sh> 
```

### automated pipeline execution starting with fastq files
```bash
# make replication timing profiles from early and late fastq files using 5000-bp window sizes and 12 threads
docker run -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq make_rt.sh \
  genome.fa 5000 12 sample1_early.fastq,sample2_early.fastq sample1_late.fastq,sample2_late.fastq
```

### step-by-step workflow

#### initial QC
```
fastqc.sh INPUT
```
#### clip and reqc
```
trim_adapters.sh INPUT
```
#### align, sort, get stats
```
align_single_end.sh INDEX NTHREADS INPUT
```
#### remove duplicates and get stats
```
remove_duplicates.sh INPUT
```
#### make windows
```
make_windows.sh CHROMSIZESFILE WINDOWSIZE
```
#### calcaulte coverage
```
window_density.sh WINDOWSFILE INPUT
```
#### filter windows
```
# filter a defined set of bedgraph files
filter_windows.R SAMPLE1_E.bg SAMPLE1_L.bg SAMPLE2_E.bg SAMPLE2_L.bg SAMPLEn_E.bg SAMPLEn_L.bg

# or use a wildcard
filter_windows.R *_[EL]_*_w5000.bg
```
#### calculate log2 ratios
```
calc_log2ratio.sh EARLYBG LATEBG
```
#### create a reference distribution from multiple files
```
make_average_distribution.sh PREFIX BG1 BG2 BGn
```
#### quantile normalize
```
quantile_normalize.sh INPUT REFERENCE
```
#### loess smooth profiles
```
loess_smooth.R SPANSIZE BG1 BG2 BGn
```













