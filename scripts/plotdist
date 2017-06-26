#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
bgfiles <- args

removeext <- function( filenames ){
	filenames<-as.character(filenames)
	for(i in 1:length(filenames)){
		namevector<-unlist(strsplit(filenames[i],"\\."))
		filenames[i]<-paste(namevector[1:(length(namevector)-1)],collapse=".")
	}
	filenames
}

for(i in seq_len(num(args)){
  outname <- paste0(removeext(basename(i)),"_hist.pdf")
  bg <- read.table(i,sep="\t",stringsAsFactors=F,header=F)
  pdf(outname)
  plot(density(bg[,4]))
  dev.off()
}
