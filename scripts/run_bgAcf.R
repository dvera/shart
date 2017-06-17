#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
df = read.table(args[1], stringsAsFactors=FALSE)
cat(acf(df[,4],lag.max=1,plot=F)[[1]][2],"\n")
