#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

EARLYFASTQ=$1
LATEFASTQ=$2
GENOMEFA=$3
WINSIZES=$4
NTHREADS=$5

CHROMSIZES=$(basename $GENOMEFA).chrom.sizes
GENOMEPREFIX=$(basename $GENOMEFA | sed 's/\.fa$//g' | sed 's/\.fasta$//g')
echo "starting repli-seq pipeline..."
echo ""
echo "early fastq is $EARLYFASTQ"
echo "late fastq is $LATEFASTQ"
echo "genome fasta is $GENOMEFA"
echo "window size is $WINSIZES"
echo "number of threads is $NTHREADS"
echo ""
echo "cancel now if what is described above is undesirable..."
sleep 10

if [ ! -f ${GENOMEPREFIX}_idx.tar.gz ]; then
  echo "creating bwa index"
  make_index.sh $GENOMEFA $GENOMEPREFIX
else
  echo "bwa index exists"
fi

echo "making chromsizes file"
make_chromsizes.sh $GENOMEFA $GENOMEPREFIX

EPREFIX=$(basename $EARLYFASTQ | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')
LPREFIX=$(basename $LATEFASTQ | sed 's/\.gz$//g' | sed 's/\.fq$//g' | sed 's/\.fastq$//g')

echo "clipping adapter sequences from reads"
${DIR}/run_cutadapt.sh $EARLYFASTQ $EPREFIX
${DIR}/run_cutadapt.sh $LATEFASTQ $LPREFIX

echo "aligning reads to genome"
${DIR}/run_bwa-mem-single.sh ${EPREFIX}_clip.fastq ${GENOMEPREFIX}_idx.tar.gz ${EPREFIX} $NTHREADS
${DIR}/run_bwa-mem-single.sh ${LPREFIX}_clip.fastq ${GENOMEPREFIX}_idx.tar.gz ${LPREFIX} $NTHREADS

echo "removing pcr duplicates"
${DIR}/run_samtools-rmdup.sh ${EPREFIX}.bam ${EPREFIX}
${DIR}/run_samtools-rmdup.sh ${LPREFIX}.bam ${LPREFIX}

echo "making windows from genome"
${DIR}/run_bedtools-makewindows.sh ${GENOMEPREFIX}.chom.sizes $WINDOWSIZE ${GENOMEPREFIX}

${DIR}/run_bedtools-intersect.sh ${EPREFIX}_rmdup.bam ${GENOMEPREFIX}_w${WINDOWSIZE}.bed ${EPREFIX}
${DIR}/run_bedtools-intersect.sh ${LPREFIX}_rmdup.bam ${GENOMEPREFIX}_w${WINDOWSIZE}.bed ${LPREFIX}

${DIR}/run_log2ratio.sh ${EPREFIX}.bg ${LPREFIX}.bg ${EPREFIX}_log2ratio





