\name{var.test.homogeneity}
\alias{var.test.homogeneity}
\title{function to perfome variance homogeneity test of different genotypic groups}
\description{
	
}
\usage{
var.test.homogeneity(trait="", formula="", trait_df, gwa_data, dir, analysis_name_specific, how_many_sigma_drop=NULL, top_snp_num_figure=3)
}
\arguments{
  \item{trait}{Name of trait you are analysing (data.frame trait_df has to contain a column with this name). If it absenses or equals to "" than the trait name must be specifed in "formula"}
  \item{formula}{Object of a class "formula". It is used in case of any adljustment nessacarity. If it absenses or equals to "" than trait must be specifed}
  \item{trait_df}{Object of class "data.frame". It contains columnes with analysing trait and covariates. Should contain column with id names.}
  \item{gwa_data}{Object of class "gwaa.data"}
  \item{dir}{Directory name where all results will be put. It should be created befored analysis starts.}
  \item{analysis_name_specific}{Any set of characters. It figures as a part of a file names in output.}
  \item{how_many_sigma_drop}{Cut off point to exclude outliers. Normally distributed trait is assumed.}
  \item{top_snp_num_figure}{For how many SNPs to make boxplot.}
}
\details{
var.test.homogeneity function perform analisis of variance homogeneity for easch SNP in data "gwa_data". The main result is pval of variance homogeneity test for each SNP.
Bartlett's test is currently used. This assumes normality in a trait distribution. Even a couple of outliers can spoil a whole pictures. 
Input parametere "how_many_sigma_drop" should be used to protect overestimation of pval.

}

\value{
As result the file "lambda.doc" is created which contains inflation factor, file with test statistics (name is qt1__adj__sex_age_test_results.txt 
				for example in case if you analyse "qt1" with adjustment on "sex" and "age"), file those contains id names which are excluded due to
	 	normality resoans, manhettan plot, two qqplots (for chi2 amd pval) and boxplots for top SNPs.
	}
  

%\references{ ~put references to the literature/web site here ~ }
\author{Maksim Struchalin}
%\note{ ~~further notes~~ 
%}
\seealso{
	\code{\link{var.meta}},
}
\examples{
	data(srdta)
	mytrait_df <- srdta@phdata[,c("id","qt1","sex","age")]
	mydir <- "mydir_var_analysis"
	dir.create(mydir)
	var.test.homogeneity(formula=qt1 ~ sex + age, 
					trait_df=mytrait_df, gwa_data=srdta, dir=mydir, analysis_name_specific="test", how_many_sigma_drop=3)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
