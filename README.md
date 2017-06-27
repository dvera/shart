# docker-4dn-hic

This repo contains the source files for a docker image stored in the docker hub container vera/docker-4dn-repliseq

## what

This repository contains a dockerfile and scripts in order to execute generate replication timing profiles from a set of raw reads from sequencing of either early- and late-replicating DNA, or from DNA extracted from cells sorted for S or G1 DNA content.

Sample data files that can be used for testing the tools are included in the `sample_data` folder.

The scripts for executing the pipeline are under the `scripts` directory and follow naming conventions `run_xx.sh`. These wrappers are copied to the docker image at build time and may be used as a single step in a workflow.

A docker image for executing these scripts can be built yourself or pulled from docker hub (vera/docker-4dn-repliseq). Images built with the dockerfile will contain both the scripts and sample data for running/testing the pipeline.

## how

### example usage
```bash
# execute a step on data in the current directory
docker run -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq <name_of_script> <args> 
```

### automated pipeline execution starting with fastq files
```bash
# make smoothed and normalized replication timing profiles from early and late fastq files using 5000-bp window sizes and 12 threads
docker run -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq repliseq  \
  genome.fa 5000 12 sample1_early.fastq,sample2_early.fastq sample1_late.fastq,sample2_late.fastq
```

### step-by-step workflow
```bash
# pull the pre-built image, create and enter a container inside the directory with your data
docker run --rm -it -h d4r -u $UID -w $PWD -v $PWD:$PWD:rw vera/docker-4dn-repliseq

# define number of CPU threads to use for the pipeline
NTHREADS=4

# download hg38 and make bwa index
wget -qO- http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > hg38.fa
index=$(index hg38.fa)

# define early and late fastq files, here using sample data
E=$(ls /opt/docker-4dn-repliseq/sample_data/*early*.fq.gz)
L=$(ls /opt/docker-4dn-repliseq/sample_data/*late*.fq.gz)

# clip adapters from reads
cfq=$(clip -t $NTHREADS $E $L)

# align reads to genome
bam=$(align -t $NTHREADS -i $index $cfq)
bstat=$(samstats -t $NTHREADS $bam)

# filter bams by alignment quality and sort by position
sbam=$(filtersort -t $NTHREADS $bam)
fbstat=$(samstats -t $NTHREADS $sbam)

# remove duplicate reads
rbam=$(dedup -t $NTHREADS $sbam)

# calculate RPKM bedGraphs for each set of alignments
bg=$(count -t $NTHREADS $rbam)

# filter windows with a low average RPKM
fbg=$(filter -t $NTHREADS $bg)

# define early and late bedGraphs
ebg=$(echo $fbg | tr ' ' '\n' | head -n $(echo $E | wc -w) | tr '\n' ',' | sed 's/,$//g')
lbg=$(echo $fbg | tr ' ' '\n' | tail -n $(echo $E | wc -w) | tr '\n' ',' | sed 's/,$//g')

# calculate log2 ratios between early and late
l2r=$(log2ratio $ebg $lbg)

# quantile-normalize replication timing profiles to the example reference bedGraph
l2rn=$(normalize -r /opt/docker-4dn-repliseq/sample_data/reference.bg $l2r)

# loess-smooth profiles using a 300kb span size
l2rs=$(smooth 300000 $NTHREADS $l2rn)
```

