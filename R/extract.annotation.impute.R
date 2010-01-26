#' extracts IMPUTE annotation
#' 
#' extracts SNP annotation information from IMPUTE files
#' 
#' 

extract.annotation.impute <- function(genofile,infofile,chromosome=NA,
		order_geno_snp_a0_a1=c(2,4:5),skip_geno=0,
		order_info_snp_pos_freq1_info_qual_type=c(2:7),skip_info=1
)
{
	if (!file.exists(genofile)) stop(paste("genofile",genofile,"does not exist"))
	if (!file.exists(infofile)) stop(paste("infofile",infofile,"does not exist"))
	
	# rsname position Al1 Al2
	if (skip_geno>0)
		tmp <- extract_text_file_columns(genofile,order_geno_snp_a0_a1)[-c(1:skip_geno),]
	else
		tmp <- extract_text_file_columns(genofile,order_geno_snp_a0_a1)
	genoannot <- as.data.frame(tmp,stringsAsFactors=FALSE)
	names(genoannot) <- c("name","A0","A1")
	
	#1snp_id 2rs_id 3position 4exp_freq_a1 5info 6certainty 7type info_type1 concord_type1 r2_type1 info_type0 concord_type0 r2_type0
	#rs_id position exp_freq_a1 info certainty type
	if (skip_info>0)
		tmp <- extract_text_file_columns(infofile,order_info_snp_pos_freq1_info_qual_type)[-c(1:skip_info),]
	else 
		tmp <- extract_text_file_columns(infofile,order_info_snp_pos_freq1_info_qual_type)
	
	infoannot <- as.data.frame(tmp,stringsAsFactors=FALSE)
	rm(tmp);gc()
	names(infoannot) <- c("name","pos","Freq1","Rsq","Quality","type")
	#print(infoannot[1:5,])
	# merge annotation
	if (!all(genoannot$name %in% infoannot$name)) {
		warning("not all snps are in info-file")
		warning(paste("missing snps:",genoannot$name[which(!(genoannot$name %in% infoannot$name))]))
	}
	mrg <- merge(genoannot,infoannot,by="name",all.x=TRUE,all.y=F)
	rownames(mrg) <- mrg$name
	mrg <- mrg[as.character(genoannot$name),]
	class(mrg$pos) <- class(mrg$Freq1) <- class(mrg$Rsq) <- class(mrg$Quality) <- "numeric"
	
	return(mrg)
}