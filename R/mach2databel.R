#'
#' converts MACH-imputed files to DatABEL (filevector) format
#' 
#' This function converts mach-imputed files to \code{DatABEL} (filevector) format.
#' After conversion, two files (outfile.fvi and outfile.fvd), corresponding 
#' to single filevector object, will appear on the disk; databel_filtered_R 
#' object connected to these files will be returned to R
#' 
#' @param mldosefile MACH mldose file name
#' @param mlinfofile MACH mlinfo file name
#' @param outfile output file name
#' 
#' @return databel_filtered_R-class object
#' 
#' @author Yurii Aulchenko
#' 
#' @keywords IO manip
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
		#print(tmp[1:10])
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
	#print(subjs[1:10])
	subjs <- sub("[0-9]*->","",subjs)
	#print(subjs[1:10])
	#print(dim(dfaobj))
	#print(length(subjs))
	dimnames(dfaobj) <- list(subjs,dnames[[2]])
	#print(dimnames(dfaobj)[[1]][1:5])
	
	if (tmpname != "") unlink(paste(tmpname,"*",sep=""))
	
	disconnect(dfaobj)
	connect(dfaobj)
	return(dfaobj)
}