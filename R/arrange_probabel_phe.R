arrange_probabel_phe <- function(modelterms,phedata,gendata,file="probabel.PHE")
{
	if (class(phedata) == "gwaa.phedata") phedata <- phphedata(phedata)
	else if (is.character(phedata)) 
		phedata <- read.table(phedata,head=T,strings=F)	
	else if (class(phedata) != "data.frame" && class(phedata) != "matrix")
		stop("phedata should be of class 'gwaa.phedata', 'data.frame', or 'matrix'")
	
	if (missing(modelterms)) 
		stop("you need to specify model terms (vector of names)")
	
	if (!any(colnames(phedata)=="id")) 
		stop("phedata should contain 'id'")
	
	if (any(is.na(match(modelterms,colnames(phedata))))) 
		stop("some modelterms are missing from phedata")
	
	if (is.character(gendata)) gendata <- databel_filtered_R(gendata)
	else if (class(gendata) == "databel_base_R") gendata <- databel_filtered_R(gendata)
	else if (class(gendata) == "databel_filtered_R") {}
	else stop("gendata should be 'databel_base/filtered_R' or name of file")
	
	dosephe <- data.frame(id=dimnames(gendata)[[1]],stringsAsFactors=F)
	phe <- merge(dosephe,phedata,by="id",all.x=T,all.y=F)
	rownames(phe) <- phe$id
	phe <- phe[as.character(dosephe$id),modelterms]
	write.table(phe,file=file,quote=F,row.names=F,col.names=T)
	
}