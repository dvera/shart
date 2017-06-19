# docker-4dn-hic

This repo contains the source files for a docker image stored in the docker hub container vera/docker-4dn-repliseq

## what

This repository contains a dockerfile and scripts in order to execute generate replication timing profiles from a set of raw reads from sequencing of either early- and late-replicating DNA, or from DNA extracted from cells sorted in S or G1 phase.

Sample data files that can be used for testing the tools are included in the `sample_data` folder. These data are not included in the docker image.

The scripts for executing the pipeline are under the `scripts` directory and follow naming conventions `run_xx.sh`. These wrappers are copied to the docker image at build time and may be used as a single step in a workflow.

The docker container for executing these scripts can be built yourself or pulled from docker hub (vera/docker-4dn-repliseq).

## how

### example usage
```bash
# execute a step on data in the current directory
docker run -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq <name_of_script> <args> 
```

### automated pipeline execution starting with fastq files
```bash
# make replication timing profiles from early and late fastq files using 5000-bp window sizes and 12 threads
docker run -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq make_rt.sh \
  genome.fa 5000 12 sample1_early.fastq,sample2_early.fastq sample1_late.fastq,sample2_late.fastq
```

### step-by-step workflow
```bash
# create and enter a container inside the directory with your data
docker run -it -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq

# initial QC
fastqc.sh INPUT

# clip and reqc
trim_adapters.sh INPUT

# align, sort, get stats
align_single_end.sh INDEX NTHREADS INPUT

# remove duplicates and get stats
remove_duplicates.sh INPUT

# make windows
make_windows.sh CHROMSIZESFILE WINDOWSIZE

# calcaulte coverage
window_density.sh WINDOWSFILE INPUT

# filter a defined set of bedgraph files
filter_windows.R *_[EL]_*_w5000.bg

# calculate log2 ratios
calc_log2ratio.sh EARLYBG LATEBG

# create a reference distribution from multiple files
make_average_distribution.sh PREFIX BG1 BG2 BGn

# quantile normalize
quantile_normalize.sh INPUT REFERENCE

# loess smooth profiles
loess_smooth.R SPANSIZE BG1 BG2 BGn
```













