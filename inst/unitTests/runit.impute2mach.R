### --- Test setup ---

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library("RUnit")
	library("GenABEL")
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

source("../inst/unitTests/shared_functions.R")

### --- Test functions ---

test.impute2mach <- function()
{
	
	library(GenABEL)
	unlink("tmp*.fv?")
	
	makedose <- function(prob) {
		dose <- 2*prob[c(F,F,T)]+prob[c(F,T,F)]
		bp <- prob[c(T,F,F)]
		miss <- which(abs(bp)<1e-16 & abs(dose)<1e-16)
		if (length(miss) > 0 ) dose[miss] <- NA
		return(dose)
	}
	makeprob <- function(prob) {
		dose2 <- prob[c(F,F,T)]
		dose1 <- prob[c(F,T,F)]
		out <- matrix(c(dose2,dose1),ncol=2)
		out
	}
	
	
	tmp0 <- read.table("tmpTEST10x15.geno",head=F)
	tmp0 <- tmp0[,c(6:dim(tmp0)[2])]
	dose <- apply(tmp0,FUN="makedose",MAR=1)
	prb <- makeprob(tmp0[1,])
	for (i in 2:dim(tmp0)[1])
		prb <- cbind(prb,makeprob(tmp0[i,]))
	dm <- dim(prb)
	prb <- as.numeric(prb)
	prb <- matrix(prb,ncol=dm[2])
	
	tmp1 <- impute2databel(geno="tmpTEST10x15.geno",
			sample="impute.sample",
			out="tmpTEST10x15_T.geno",
			old=TRUE)
	tmp1
	tmp1_m <- as(tmp1,"matrix")
	tmp1_m
	
	
	tmp0[1:5,1:6]
	tmp2 <- impute2databel(geno="tmpTEST10x15.geno",
			sample="impute.sample",
			out="tmpTEST10x15_F.geno",
			old=FALSE)
	tmp2
	tmp2_m <- as(tmp2,"matrix")
	
	
	tmp3 <- databel_filtered_R("tmpTEST10x15_F.geno.prob")
	tmp3_m <- as(tmp3,"matrix")
	prb
	tmp3_m
	
	table(abs(dose-tmp1_m)<1e-6)
	table(abs(dose-tmp2_m)<1e-6)
	table(abs(prb-tmp3_m)<1e-6)
	table(abs(tmp1_m-tmp2_m)<1e-8)
	identical(tmp1_m,tmp2_m)
	
	checkEqualsNumeric(dose,tmp1_m,tolerance=4*sqrt(.Machine$double.eps))
	checkEqualsNumeric(dose,tmp2_m,tolerance=4*sqrt(.Machine$double.eps))
	checkEqualsNumeric(prb,tmp3_m,tolerance=4*sqrt(.Machine$double.eps))
	checkIdentical(tmp1_m,tmp2_m)
	
	unlink("tmp*.fv?")	
	
}

