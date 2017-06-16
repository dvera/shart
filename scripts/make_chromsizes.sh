#!/usr/bin/env bash

GENOMEFA=$1
PREFIX=$2

CURDIR=$PWD

awk '{
  if($0~">"){
    if(length(chrom)!=0){
      print chrom,clen
    };
    chrom=$1;
    clen=0
   } else{
    clen+=length($0)
   }
  } END{
    print chrom,clen
  }' OFS='\t' $GENOMEFA | tr -d '>' | sort -k1,1 > ${CURDIR}/${PREFIX}.chrom.sizes
  
  
