roxy_files <- c(
		"add.phdata.R",
		"annotation",
		"arrange_probabel_phe.R",
		"del.phdata.R",
		"export.plink.R",
		"extract.annotation.impute.R",
		"extract.annotation.mach.R",
		"GenABEL-package.R",
		"impute2databel.R",
		"impute2mach.R",
		"mach2databel.R",
		"phdata.R",
		"polygenic.R",
		"summary.scan.gwaa.R"
)

library(roxygen)
setwd("R")
unlink("GenABEL",recursive=TRUE)
package.skeleton("GenABEL",code_files=roxy_files)
roxygenize("GenABEL",roxygen.dir="GenABEL",copy=F,unlink=F)
lstf <- paste("GenABEL/man/",list.files("GenABEL/man/"),sep="")
lstf
file.copy(lstf,"../man/",overwrite=TRUE)
unlink("GenABEL",recursive=TRUE)
