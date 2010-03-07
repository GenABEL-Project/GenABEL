#'
#' converts IMPUTE-imputed files to DatABEL (filevector) format
#' 
#' this function converts IMPUTE-imputed files to DatABEL (filevector) format
#' containing estimated dosages. 
#' After conversion, two files (outfile.fvi and outfile.fvd), corresponding 
#' to single filevector object, will appear on the disk; databel_filtered_R 
#' object connected to these files will be returned to R.
#' 
#' @param genofile IMPUTE genotype file name
#' @param samplefile IMPUTE sample file name
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
	tmpname <- get_temporary_file_name()
	tmp_fv <- text2filevector(infile=genofile,outfile=tmpname,
			rownames=2,skipcols=5,
			#skiprows,
			transpose=TRUE,R_matrix=FALSE,type="FLOAT")
	if ((dim(tmp_fv)[1] %% 3) != 0) 
		stop("strange data in IMPUTE geno file: number of columns - 5 not dividable by 3")
	
	makedose <- function(prob) {
		dose <- 2*prob[c(F,F,T)]+prob[c(F,T,F)]
		bp <- prob[c(T,F,F)]
		miss <- which(abs(bp)<1e-16 & abs(dose)<1e-16)
		if (length(miss) > 0 ) dose[miss] <- NA
		return(dose)
	}
	
	pfun <- function(a) return(a)
	
	#dosefile <- make_empty_fvf(paste(outfile,".dose",sep=""), 
	#		nvariables = dim(tmp_fv)[2], nobservations = round(dim(tmp_fv)[1]/3),type = "FLOAT")
	
	#print(dimnames(tmp_fv)[[2]])
	#print(dim(tmp_fv))
	#print(class(tmp_fv))
	saved_names <- get_dimnames(tmp_fv)
	#print(class(tmp_fv))
	
	tmplst <- list( as.character(1:dim(tmp_fv)[1]) , as.character(1:dim(tmp_fv)[2]) )
	#print(class(tmp_fv))
	#print(dimnames(tmp_fv)[[2]])
	dimnames(tmp_fv) <- tmplst
	#print(class(tmp_fv))
	#print(dimnames(tmp_fv)[[2]])
	
	#print("before apply2dfo")
	#print("calling apply2dfo")
	dosefile <- apply2dfo(dfodata=tmp_fv, anFUN = "makedose", 
			MAR = 2, procFUN = "pfun",prob=SNP,
			outclass="databel_filtered_R",
			outfile=paste(outfile,".dose",sep=""),
			type="FLOAT",transpose=FALSE)
	#print("after apply2dfo")
	dimnames(dosefile)[[2]] <- saved_names[[2]]
	#print("dimnames [[2]]")
	
	if (!missing(samplefile))
	{
		temp <- scan(samplefile,what="character",nlines=1)
		samnames <- scan(samplefile,what="character",skip=2)
		samnames <- samnames[c(F,T,rep(F,(length(temp)-2)))]
		if (length(samnames) == dim(dosefile)[1]) 
			rownames(dosefile) <- samnames
		else
			warning("number of IDs specified in sample file does not match to geno-dimension; dropping ID names")
	} else 
		warning("sample file not specified, you will not be able to use ID names (only index)")
	
	rm(tmp_fv);gc();unlink(paste(tmpname,"*",sep=""))
	
	disconnect(dosefile)
	connect(dosefile)
	return(dosefile)
}