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

test.convert.snp <- function()
{
	library(GenABEL)
	convert.snp.tped(tped="gatest.tped", tfam="gatest.tfam", out="gatest.raw", strand="+")
	data <- load.gwaa.data(gen="gatest.raw", phe="gatest.phe")
	data
	coding(data)
	as.character(data)
	as.numeric(data)
	convert.snp.tped(tped="gatest1.tped", tfam="gatest1.tfam", out="gatest1.raw", strand="+")
	data1 <- load.gwaa.data(gen="gatest1.raw", phe="gatest1.phe")
	data1
	coding(data1)
	as.character(data1)
	as.numeric(data1)
	mrg <- merge(gtdata(data),gtdata(data1))
	mrg
	coding(mrg$data)
	as.character(mrg$data)
	as.numeric(mrg$data)
	
	coding(data)
	class(data)
	coding(data) <- coding(data1)
	class(data)
	coding(data)
	coding(data1)
	as.character(data)
	as.character(data1)
	mrg <- merge(gtdata(data),gtdata(data1))
	mrg
	coding(mrg$data)
	as.character(mrg$data)
	as.numeric(mrg$data)
	unlink("*.raw")
}