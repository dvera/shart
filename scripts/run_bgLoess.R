#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
lspan <- as.numeric(args[1])
bgfiles <- args[2:length(args)]

removeext <- function( filenames ){
	filenames<-as.character(filenames)
	for(i in 1:length(filenames)){
		namevector<-unlist(strsplit(filenames[i],"\\."))
		filenames[i]<-paste(namevector[1:(length(namevector)-1)],collapse=".")
	}
	filenames
}

	options(scipen=9999)
	# assumes sorted bg
	if(lspan<1){
		stop("lspan must indicate distance in bp")
	}
	bgnames <- basename(removeext(bgfiles))
	numbgs <- length(bgfiles)
	outnames <- paste(bgnames,"_loess",lspan,".bg",sep="")
	
	#dump <- mclapply(seq_len(numbgs), function(x){
	dump <- lapply(seq_len(numbgs), function(x){
		curbg <- as.data.frame( read.table ( bgfiles[x], header=FALSE, stringsAsFactors=FALSE ) )
		
		chroms    <- unique(curbg[,1])
		numchroms <- length(chroms)
		windowsize <- median(curbg[,3]-curbg[,2])
		all=split(curbg,curbg[,1])
		
		#lscores<-mclapply(1:numchroms,function(i){
		lscores<-lapply(1:numchroms,function(i){
			
			cur <- all[[i]]
			chromlspan <- lspan/sum(cur[,3]-cur[,2])
			
			cura <- as.data.frame(lapply(4:ncol(cur), function(k){
				cur[,k] <- tryCatch({
						loess(cur[,k]~cur[,2],span=chromlspan)$fitted
					},warning = function(war){
						print(paste("warning for file",bgnames[x],"chromosome",cur[1,1],":",war))
						out <- loess(cur[,k]~cur[,2],span=chromlspan)$fitted
						return(out)
					},error = function(err){
						print(paste("smoothing failed for file",bgnames[x],"chromosome",cur[1,1],":",err))
						return(cur[,k])
					}
				)
			}))
			
			cura <- cbind(cur[,1:3],cura)
			colnames(cura) <- colnames(cur)
			return(cura)
		})
		smoothbg<-do.call(rbind,lscores)
		write.table(smoothbg,file=outnames[x],col.names=FALSE,quote=FALSE, sep="\t", row.names=F)
		rm(curbg)
		gc()
	})
	
