#!/usr/bin/env bash

GENOMEFA=$1
PREFIX=$2

CURDIR=$PWD

bwa index -p $PREFIX $GENOMEFA
FADIR=$(readlink -f $GENOMEFA)

IDXFILES="${PREFIX}.amb ${PREFIX}.ann ${PREFIX}.bwt ${PREFIX}.pac ${PREFIX}.sa"
(cd $FADIR && tar -czf $CURDIR/${PREFIX}_idx.tar.gz $IDXFILES && rm $IDXFILES)
