#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
df = read.table(args[1], stringsAsFactors=FALSE)
write(acf(df[,4],lag.max=1,plot=F)[[1]][2],paste0(args[1],".acf"))

