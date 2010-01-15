#=====================================================================================
#
#       Filename:  homo_var_analysis.R 
#
#    Description:  Analizis of homogeneity variance
#
#        Version:  1.1
#        Created:  14-Apr-2009
#       Revision:  none
#				Modified:  22-Apr-2009
#       
#
#         Author:  Maksim V. Struchalin
#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
#          Email:  m.struchalin@erasmusmc.nl
#
#=====================================================================================


"var.test" <-
function(trait="", formula="", trait_df, gwa_data, dir, analysis_name_specific) {




if(class(trait) == "formula")
	{
	stop("trait has to have type character")
	}


if(trait=="" && formula=="")
	{
	stop("You have bot specifed trait or residuals which has to be analysed")
	}


if(formula != "" && class(formula) != "formula")
	{
	message <- paste("formula has invalid type (",class(formula), ")", sep="" )
	stop(message)
	}


if(class(trait_df) != "data.frame")
	{
	message <- paste("trait_df has invalid type (", class(trait_df), ")", sep="")
	stop(message)
	}


if(class(gwa_data) != "gwaa.data")
	{
	message <- paste("gwa_data has invalid type (", class(gwa_data), ")", sep="")
	stop(message)
	}

if(length(unique(colnames(trait_df) == "id")) == 1)
	{
	message <- paste("gwa_data doesn't have column name id. Only ", colnames(trait_df))
	stop(message)
	}


analysis_name <- "homo_variance"

if(formula != "")
	{
	formula_terms <- as.character(terms(formula))
	traits <- c(formula_terms[2], strsplit(formula_terms[3], " \\+ ", perl=T)[[1]])
	traits_num <- length(traits)
	analysis_name <- paste(traits[1], "__adj_", sep="")
	for(i in 2:traits_num)
		{
		analysis_name <- paste(analysis_name, traits[i], sep="_")
		}


	trait <- traits[1]
	}
else
	{
	analysis_name <- trait
	}

analysis_name <- paste(analysis_name, analysis_name_specific, sep="_")





#library(GenABEL)
		
#source("/home/maksim/work/GenABEL_dev/GenABEL/R/var.homogeneity.test.R")
#dyn.load("/home/maksim/work/GenABEL_dev/GenABEL/src/var.homogeneity.test.so")

#source("~/work/mylib/gwa.plot.R")
#source("~/work/mylib/qqplot.R")

#source("/home/maksim/work/mylib/drop.out.points.R")
#source("/home/maksim/work/mylib/get_snp_plots.R")



trait_names <- c(trait, strsplit(as.character(terms(formula))[3], " \\+ ", perl=T)[[1]])
trait_df_nona_ <- na.omit(trait_df[, c("id", trait_names)])

#let's drop out points which are beyond 3*sigma
trait_df_nona_[,trait] <- drop.out.points(trait_df_nona_[,trait])


excluded_ids <- trait_df_nona_$id[is.na(trait_df_nona_[,trait])]
#write.table(data.frame(ids=excluded_ids), file=paste(dir, "/", analysis_name, "_excluded_ids.txt", sep=""), row.names=F, quote=F)

print("excluded ids bacause of normality reason:")
print(excluded_ids)

trait_df_nona <- na.omit(trait_df_nona_)



trait_df_nona$residuals <- lm(formula, data=trait_df_nona)$residuals





bitmap(paste(dir, "/", analysis_name, "_", analysis_name_specific, "resid_hist.jpeg", sep=""), type="jpeg")
hist(trait_df_nona$residuals, main=analysis_name)
dev.off()

data <- add.phdata(gwa_data, trait_df_nona[,c("id", "residuals")])
print(paste("We have", length(trait_df_nona$residuals), "ids"))



var_homo_chi2_trait_residuals <- var.homogeneity.test(data@phdata$residuals, data)



var_homo_chi2_trait_residuals_df <- data.frame(snpname=as.character(data@gtdata@snpnames) ,pval=pchisq(var_homo_chi2_trait_residuals, df=2, lower.tail=F), 
				position=data@gtdata@map, chromosome = data@gtdata@chromosome)


var_homo_chi2_trait_residuals_df_sorted <- var_homo_chi2_trait_residuals_df[order(var_homo_chi2_trait_residuals_df$pval),]
write.table(var_homo_chi2_trait_residuals_df_sorted, file=paste(dir, "/", analysis_name, "_results.txt", sep=""), row.names=F, quote=F)

#print(dim(var_homo_chi2_trait_residuals_df))
#print(var_homo_chi2_trait_residuals_df[1:5,])
print(paste(dir, "/", analysis_name, "_var_manhettan", sep=""))

gwa.plot(var_homo_chi2_trait_residuals_df, filename=paste(dir, "/", analysis_name, "_var_manhettan", sep=""))
qqplot.plot(var_homo_chi2_trait_residuals, df=2, filename=paste(dir, "/", analysis_name, "_var_homo", sep=""))


lambda <- median(var_homo_chi2_trait_residuals, na.rm=T)/qchisq(0.5, df=2)

system(paste("echo lambda=", lambda, " > ", dir, "/lambda.doc", sep=""))
print(paste("lambda=", lambda, sep=""))


#plot top 20 SNPs


gen <- data.frame(id=as.character(data@gtdata@idnames))
rownames(gen) <- gen$id

for(snp in 1:50)
	{
	snpname <- as.character(var_homo_chi2_trait_residuals_df_sorted$snpname[snp])
	print(paste(snp , ": ", snpname, sep=""))
	gen[,snpname] <- as.numeric(data[,snpname])
#	bitmap(paste(dir, "/", analysis_name, "-", snpname, ".jpeg", sep=""), type="jpeg")
	
	gen_trait_df <- data.frame(id=trait_df_nona$id ,trait=trait_df_nona$residuals)

	gen_trait_df[,snpname] <- as.numeric(gen[as.character(gen_trait_df$id),snpname])
	
	gen_trait_df_nona <- na.omit(gen_trait_df)
	
	chr <- as.character(var_homo_chi2_trait_residuals_df_sorted$chromosome[snp])
	pval <- var_homo_chi2_trait_residuals_df_sorted$pval[snp]

	get.snp.distr.plot(df=data.frame(trait=gen_trait_df_nona$trait, snp=gen_trait_df_nona[,snpname]), plotname=analysis_name, snpname=snpname, chromosome=chr,
				 						 traitname=trait, filename=paste(dir, "/", analysis_name, "_",snpname, ".jpeg", sep=""))

#	plot(gen_trait_df_nona[,snpname], gen_trait_df_nona$trait, main=paste(analysis_name, ",\n", snpname, ", chr=", chr,
#			  ", pval=", pval,  sep=""), xlab=snpname, ylab=analysis_name)
#	dev.off()
	}





save.image(paste(dir, "/",analysis_name,"_image.RData", sep=""))





}
