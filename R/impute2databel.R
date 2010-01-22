#'
#' converts IMPUTE-imputed files to DatABEL (filevector) format
#' 
#' this function converts IMPUTE-imputed files to DatABEL (filevector) format
#' computing dosages
#' 
#' @param genofile MACH mldose file name
#' @param samplefile MACH mlinfo file name
#' @param outfile output file name
#' 
#' @return databel_filtered_R-class object
#' 
#' 


impute2databel <- function(genofile,samplefile,outfile) # dosefile = TRUE 
{
	if (!require(DatABEL))
		stop("this function requires DatABEL package to be installed")
	if (missing(genofile))
		stop("genofile file must be specified")
	if (missing(outfile)) outfile <- genofile
# extract snp names (varnames)
	#   if (tmpname != "")
	#       text2filevector(infile=genofile,outfile=outfile,
	#               colnames=tmpname,
	#               rownames=2,skipcols=5,
	#               #skiprows,
	#               transpose=TRUE,R_matrix=FALSE,type="FLOAT")
	#   else 
	tmpname <- paste("tmp",round(runif(1,min=10000,max=1000000)),sep="")
	while (file.exists(tmpname)) 
	{
		tmpname <- paste("tmp",round(runif(1,min=10000,max=1000000)),sep="")
	}
	tmp_fv <- text2filevector(infile=genofile,outfile=tmpname,
			rownames=2,skipcols=5,
			#skiprows,
			transpose=TRUE,R_matrix=FALSE,type="FLOAT")
	if ((dim(tmp_fv)[1] %% 3) != 0) 
		stop("strange data in IMPUTE geno file: number of columns - 5 not dividable by 3")
	
	makedose <- function(prob) return(2*prob[c(F,F,T)]+prob[c(F,T,F)])
	pfun <- function(a) return(a)
	
	#dosefile <- make_empty_fvf(paste(outfile,".dose",sep=""), 
	#		nvariables = dim(tmp_fv)[2], nobservations = round(dim(tmp_fv)[1]/3),type = "FLOAT")
	
	dosefile <- apply2dfo(dfodata=tmp_fv, anFUN = "makedose", 
			MAR = 2, procFUN = "pfun",prob=SNP,
			outclass="databel_filtered_R",
			outfile=paste(outfile,".dose",sep=""),
			type="FLOAT",transpose=FALSE)
	
	
	if (!missing(samplefile))
	{
		samnames <- scan(samplefile,what="character",skip=1)
		rownames(dosefile) <- samnames
	} else 
		warning("sample file not specified, you will not be able to use ID names (only index)")
	
	rm(tmp_fv);gc();unlink(tmpname)
	
	return(dosefile)
}