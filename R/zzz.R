.onLoad <- function(lib, pkg) {
    GenABEL.version <- "1.6-5"
    cat("GenABEL v.",GenABEL.version,"(September 4, 2010) loaded\n")
    
    address <- c(
            "http://genabel.r-forge.r-project.org/version_and_news.html",
            "http://mga.bionet.nsc.ru/~yurii/ABEL/version_and_news.html"
    )
    svtmo <- options("timeout")
    options("timeout"=10)
    tryRes1 <- 0; class(tryRes1) <- "try-error"
    curaddr <- 1
    while (class(tryRes1) == "try-error" && curaddr <= length(address) ) {
        suppressWarnings(
                tryRes0 <- try(conn <- url(address[curaddr]),silent=TRUE)
        )
        suppressWarnings(
                tryRes1 <- try(fulltext <- readLines(conn),silent=TRUE)
        )
		close(conn)
		curaddr <- curaddr + 1
    }
    if (class(tryRes1) != "try-error") {
        if (length(fulltext)>0)
        {
			# compare versions
            a <- tolower(fulltext)
            a <- a[grep("<gastable>",a)+1]
            if (length(a)>0) {
                a <- strsplit(a,"")[[1]]
                ver <- a[grep("[0-9]",a)]
                ver <- paste(ver[1],".",ver[2],"-",ver[3],sep="")
                if (GenABEL.version != ver) {
                    cat(  "\nInstalled GenABEL version (",GenABEL.version,") is not the same as stable\n",
                            "version available from CRAN (",ver,"). Unless used intentionally,\n",
                            "consider updating to the latest CRAN version. For that, use\n",
							"'install.packages(\"GenABEL\")', or ask your system administrator\n",
							"to update the package.\n\n",sep="")
					strnews <- grep("<ganews>",tolower(fulltext))
					endnews <- grep("</ganews>",tolower(fulltext))
					if ((endnews-1) >= (strnews+1)) {
						cat(fulltext[(strnews+1):(endnews-1)],sep="\n")
					}
				}
            }
        }
        #rm(a,fulltext,ver)
    }
    options("timeout"=svtmo)
    #rm(tryRes0,tryRes1,conn,svtmo)
}