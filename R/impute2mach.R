#' converts IMPUTE to MACH files
#'
#' function to convert IMPUTE files to MAC format
#' 

impute2mach <- function(genofile,infofile,samplefile,machbasename, maketextdosefile = TRUE    , ... )
{
    if (!require(DatABEL))
        stop("this function requires DatABEL package to be installed")
    
    if (!is.character(machbasename)) stop("machbasename must be character")
    if (length(machbasename) == 1) {
        machdose <- paste(machbasename,".dose",sep="")
        machinfo <- paste(machbasename,".info",sep="")
        machlegend <- paste(machbasename,".legend",sep="")
    } else if (length(machbasename) == 3) {
        if (length(unique(machbasename)) != 3) stop("names must be unique")
        machdose <- machbasename[1]
        machinfo <- machbasename[2]
        machlegend <- machbasename[3]
    } else stop("machbasename must be character of length 1 or 3")
    
    # create temporary DA file
    dfo <- impute2databel(genofile=genofile,samplefile=samplefile,outfile=genofile)
    if (maketextdosefile) {
        tmpname2 <- get_temporary_file_name()    
        # transpose file
        # ... 
        dfo <- as(dfo,"matrix")
    }
    # get annotattion
    #print("AAA")
    annot <- extract.annotation.impute(genofile=genofile,infofile=infofile, ... )
    #print(annot[1:5,])
    
    # arrange MLINFO and legend file
    #SNP    Al1    Al2    Freq1    MAF    Quality    Rsq
    annot$MAF <- pmin(annot$Freq1,(1.-annot$Freq1))
    #print(annot[1:5,])
    info_annot <- annot[,c("name","A1","A0","Freq1","MAF","Quality","Rsq")]
    #print(info_annot[1:5,])
    write.table(info_annot,file=machinfo,row.names=FALSE,col.names=TRUE,quote=F,sep="\t")
    legend_annot <- annot[,c("name","pos","A1","A0")]
    #print(legend_annot[1:5,])
    write.table(legend_annot,file=machlegend,row.names=FALSE,col.names=TRUE,quote=F,sep="\t")
    
    if (maketextdosefile) {
        # arrange MLDOSE file
        ids <- dimnames(dfo)[[1]]
        if (file.exists(machdose)) unlink(machdose)
        outfile <- file(machdose,open="wt")
        for (i in 1:dim(dfo)[1])
        {
            # when using transposed DA object, use as.vector(dfo[,i]) (COLUMN!!!)
            outline <- c(ids[i],"MLDOSE",as(dfo[i,],"vector"))
            #print(outline)
            #print(i)
            #print(class(outline))
            write(x=outline,file=outfile,append=TRUE,ncolumns=length(outline),sep=" ")
            if ((i %% 100)==0 || i==dim(dfo)[1]) print(i)
        }
        close(outfile)
        unlink(paste(tmpname2,"*",sep=""))
    }
    #print("AAA")
    disconnect(dfo)
    rm(dfo);gc()
}