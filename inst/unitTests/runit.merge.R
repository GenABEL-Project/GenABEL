### --- Test setup ---
#
# regression test
#

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library(RUnit)
	library(GenABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
#source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.merge.bug1641 <- function()
{
	library(GenABEL)
	print(  packageVersion("GenABEL") )
	data(srdta)
	x1 <- srdta[1:4,1]
	x2 <- srdta[5:10, 2]
	xy <- merge(x1, x2, intersected_snps_only=FALSE)
}