#'
#' converts mach-imputed files to DatABEL (filevector) format
#' 
#' this function converts mach-imputed files to DatABEL (filevector) format
#' 
#' @param mldosefile MACH mldose file name
#' @param mlinfofile MACH mlinfo file name
#' @param outfile output file name
#' 
#' @value databel_filtered_R-class object
#' 
#' 


mach2databel <- function(mldosefile,mlinfofile,outfile) 
{
	if (!require(DatABEL))
		stop("this function requires DatABEL package to be installed")
	if (missing(mldosefile))
		stop("mldose file must be specified")
	if (missing(outfile)) outfile <- mldosefile
# extract snp names (varnames)
	tmpname <- ""
	if (!missing(mlinfofile))
	{
		tmp <- scan(mlinfofile,what="character",skip=1)
		tmp <- tmp[c(T,F,F,F,F,F,F)]
		print(tmp[1:10])
		tmpname <- get_temporary_file_name()
		write(tmp,file=tmpname)
		rm(tmp);gc()
	} else 
		warning("mlinfo file not specified, you will not be able to use snp names (only index)")
	
	if (tmpname != "")
		dfaobj <- text2filevector(infile=mldosefile,outfile=outfile,
				colnames=tmpname,
				rownames=1,skipcols=2,
				#skiprows,
				transpose=FALSE,R_matrix=FALSE,type="FLOAT")
	else 
		dfaobj <- text2filevector(infile=mldosefile,outfile=outfile,
				rownames=1,skipcols=2,
				#skiprows,
				transpose=FALSE,R_matrix=FALSE,type="FLOAT")

	dnames <- get_dimnames(dfaobj)
	subjs <- dnames[[1]]
	print(subjs[1:10])
	subjs <- sub("[0-9]*->","",subjs)
	print(subjs[1:10])
	print(dim(dfaobj))
	print(length(subjs))
	set_dimnames(dfaobj) <- list(subjs,dnames[[2]])
	print(get_dimnames(dfaobj)[[1]][1:5])
	
	if (tmpname != "") unlink(paste(tmpname,"*",sep=""))
	
	return(dfaobj)
}