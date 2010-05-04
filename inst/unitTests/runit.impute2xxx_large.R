### --- Test setup ---

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library("RUnit")
	library("GenABEL")
}
#test.empty <- function(){}
### do not run
#stop("SKIP THIS TEST")
###

# common functions and data

data_prefix = "big_data/KOM"

### --- Test functions ---

test.impute2xxx_large <- function()
{

	geno <- paste(data_prefix,".geno",sep="")
	info <- paste(data_prefix,".info",sep="")
	samp <- paste(data_prefix,".sample",sep="")
	
	if (!(file.exists(geno) && file.exists(info) && file.exists(samp))) {
		warning(paste("big files are missing; stopping tests"))
		return(NULL)
	}
	
	out <- "KOM_test_out.geno"
	fvdose <- impute2databel(geno=geno,
			sample=samp,
			out=out,
			makeprob=TRUE,
			old=FALSE)
	
}
