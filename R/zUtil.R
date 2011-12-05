# utility functions invisible to user
lossFunctionLambdaKS <- function(lambda,chi2values, ... ) {
	ksT <- ks.test(chi2values/lambda, ... )$stat
	return( ksT )
}
estLambdaKS <- function(chi2values,limits=c(0.5,100)) {
	iniLambda <- 1
	optRes <- optimize(lossFunctionLambdaKS, interval=limits, chi2values=chi2values, "pchisq", 1)
	return(optRes$minimum)
}
