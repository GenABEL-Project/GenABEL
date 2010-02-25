roxy_files <- c(
		"add.phdata.R",
		"arrange_probabel_phe.R",
		"del.phdata.R",
		"extract.annotation.impute.R",
		"extract.annotation.mach.R",
		"GenABEL-package.R",
		"impute2databel.R",
		"impute2mach.R",
		"mach2databel.R",
		"phdata.R"
)

library(roxygen)
setwd("R")
system("rm -r GenABEL")
package.skeleton("GenABEL",code_files=roxy_files)
roxygenize("GenABEL",roxygen.dir="GenABEL",copy=F,unlink=F)
system("cp GenABEL/man/*Rd ../man/.")
system("rm -r GenABEL")
