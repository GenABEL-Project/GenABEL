#' Estimation of polygenic model
#' 
#' Estimation of polygenic model using hierarchical
#' GLM (hglm package)
#' 
#' @param formula Formula describing fixed effects to be used in analysis, e.g. 
#' y ~ a + b means that outcome (y) depends on two covariates, a and b. 
#' If no covariates used in analysis, skip the right-hand side of the 
#' equation.
#' @param kinship.matrix Kinship matrix, as provided by e.g. ibs(,weight="freq"), 
#' or estimated outside of GenABEL from pedigree data.
#' @param data An (optional) object of \code{\link{gwaa.data-class}} or a data frame with 
#' outcome and covariates
#' @param family a description of the error distribution and link function to be 
#' used in the mean part of the model (see \code{\link{family}} for details of 
#' family functions)
#' @param ... other parameters passed to \code{\link{hglm}} call
#' 
#' @author Xia Shen, Yurii Aulchenko
#' 
#' @references 
#' need to put reference here
#' 
#' @seealso 
#' \code{\link{mmscore}},
#' \code{\link{grammar}}
#' \code{\link{polygenic}}
#' 
#' @examples 
#' data(ge03d2ex.clean)
#' df <- ge03d2ex.clean[,autosomal(ge03d2ex.clean)]
#' gkin <- ibs(df,w="freq")
#' 
#' # for quantitative traits
#' h2ht <- polygenic_hglm(height ~ sex + age,kin=gkin,df)
#' # estimate of heritability
#' h2ht$esth2
#' # other parameters
#' h2ht$h2an
#' 
#' # for binary traits
#' h2dm <- polygenic_hglm(dm2 ~ sex + age,kin=gkin,df,family=binomial(link = 'logit'))
#' # estimated parameters
#' h2dm$h2an
#' 
#' 
#' @keywords htest
#' 
"polygenic_hglm" <- function(formula,kinship.matrix,data,family=gaussian(), ... )
{
	if (!require(hglm))
		stop("this function requires 'hglm' package to be installed")
	if (!missing(data)) if (is(data,"gwaa.data")) 
		{
#			checkphengen(data)
			data <- phdata(data)
		}
	if (!missing(data)) 
		if (!is(data,"data.frame")) 
			stop("data should be of gwaa.data or data.frame class")
	allids <- data$id
	
	relmat <- kinship.matrix
	relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]
	mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
	y <- model.response(mf)
	desmat <- model.matrix(formula,mf)
	phids <- rownames(data)[rownames(data) %in% rownames(mf)]
	mids <- (allids %in% phids)
	relmat <- relmat[mids,mids]
	relmat <- relmat*2.0
	s <- svd(relmat)
	L <- s$u %*% diag(sqrt(s$d))
	res_hglm <- hglm(y = y, X = desmat, Z = L, family = family, ... )
	#sum_res_hglm <- summary(res_hglm)
	
	out <- list()
	
	out$measuredIDs <- mids
	out$hglm <- res_hglm
	out$h2an <- list()
	tVar <- res_hglm$varRan+res_hglm$varFix
	out$esth2 <- res_hglm$varRan/tVar
	out$h2an$estimate <- c(res_hglm$fixef,out$esth2,tVar)
	names(out$h2an$estimate)[(length(out$h2an$estimate)-1):(length(out$h2an$estimate))] <- 
			c("h2","totalVar")
	out$pgresidualY <- rep(NA,length(mids))
	out$pgresidualY[mids] <- y-res_hglm$fv
	out$residualY <- rep(NA,length(mids))
	out$residualY[mids] <- y - desmat %*% res_hglm$fixef
	out$InvSigma <- ginv(tVar*out$esth2*relmat + 
					diag(tVar*(1-out$esth2),ncol=length(y),nrow=length(y)))
	
	class(out) <- "polygenic"
	out
}
