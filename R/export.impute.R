"export.impute" <- function(data,genofile="impute.gen",samplefile="impute.sample",
		strandfile="impute.strand",cachesizeMb=128) 
{
	cat("beta-version of export.impute:\n\tstrand file not generated (what is format?)\n\tsample file contains IDs of individulas\n")
	
	if (!is(data,"gwaa.data")) 
		stop("Data argumet should be of gwaa.data-class")
	if (length(levels(data@gtdata@chromosome)) > 1) 
		warning("data contains > 1 chromosome; this can not be reflected in IMPUTE format")
	if(cachesizeMb<10) {
		warning("cachesizeMB < 10, set to 10")
		cachesizeMb <- 10
	}
	# write sample file
	samples <- matrix(data@gtdata@idnames,ncol=1)
	colnames(samples) <- "id"
	write.table(samples,file=samplefile,col.names=TRUE,quote=F)
	rm(samples)
	gc()
	
	# write genotype file
	cat("writing sample file ...")
	rsNames <- as.character(data@gtdata@snpnames)
	rsPos <- as.integer(data@gtdata@map)
	coding <- as.character(data@gtdata@coding)
	allele1=alleleID.reference()[coding]
	allele2=alleleID.effective()[coding]
	rm(coding)
	gc()
	cat("... done!\n")
	# write snp by snp, replace missing with HWE proportions
	cat("writing genotypes ...\n")
#	freqs <- summary(data@gtdata)[,"Q.2"]
	noutsnps <- ceiling((cachesizeMb*1024*1024)/(3*8*data@gtdata@nids))
#	print(noutsnps)
	bfrom <- seq(1,data@gtdata@nsnps,noutsnps)
	bto <- bfrom+noutsnps-1
	bto[length(bto)] <- data@gtdata@nsnps
#	print(bfrom)
#	print(bto)
	for (i in 1:length(bfrom)) {
		print(paste("writing SNPs from",bfrom[i],"to",bto[i],"..."))
		tmpdata <- data@gtdata[,c(bfrom[i]:bto[i])]
		tmp <- .Call("get_impute_snp_matrix",tmpdata@nids,tmpdata@nsnps,tmpdata@gtps) #,freqs)
		rm(tmpdata);gc()
		tmp <- cbind(rsNames[bfrom[i]:bto[i]],rsNames[bfrom[i]:bto[i]],rsPos[bfrom[i]:bto[i]],
				allele1[bfrom[i]:bto[i]],allele2[bfrom[i]:bto[i]],tmp)
# this is ugly: t(out)
		if (i == 1)
			write(t(tmp),file=genofile,ncolumns=(dim(tmp)[2]),append=FALSE)
		else
			write(t(tmp),file=genofile,ncolumns=(dim(tmp)[2]),append=TRUE)
	}
	cat("... done!\n")
}
