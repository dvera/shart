#!/usr/bin/env bash

CHROMSIZES=$1
WINDOWSIZE=$2
PREFIX=$3

bedtools makewindows -w $WINDOWSIZE -g $CHROMSIZES > ${PREFIX}_w${WINDOWSIZE}.bed
