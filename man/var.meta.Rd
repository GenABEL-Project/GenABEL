\name{var.meta}
\alias{var.meta}
\title{function to pool genotypic means}
\description{
	Function takes file names where trait means separatly for different genotypes are stored for each SNP and
 pool (metaanalyse) them. Input files are in format mentioned below (plink format).
}
\usage{
var.meta(filenames_array ,output_filename="output.variance.metaanalysis", exclude_snp_below_threshold=F, threshold=30, all_warnings=F)
}
\arguments{
  \item{filenames_array}{Set of files which you are going to pool}
  \item{output_filename}{Name of output file where pooled snps are stored. "output.variance.metaanalysis" by default}
  \item{exclude_snp_below_threshold}{kind of "Minor Allele Frequency Filter". But absolute is number is used (threshold)}
  \item{threshold}{If at least one genotypic group has number of id < threshold and previuos option is T than SNO is excluded}
  \item{all_warnings}{T if we interested in which snps are excluded because of filter}
}
\details{
	Function is intended for pooling (metaanalysing) trait means for different genotypes separetly. 
	Input files have to be like in example below (plink generate this kind of files).
	At the end you have one file in same format but with pooled MEANs and SDs
}

\value{
	Pooled means are stored in file with name which pointed in output\_filename variable.
	Function returns data.frame object with inflation factors (lambda) for each input file and for pooled snps.
	Inflation factors is calculated as 2df bartlett's (variance homogeneity) test. 

}
\author{Maksim Struchalin}


\examples{
\dontrun{
 first input file:
 CHR SNP        VALUE  G11    G12 G22
 1   rs12497307 GENO   A/A    A/G  G/G
 1   rs12497307 COUNTS  0      2  5299
 1   rs12497307 FREQ    0    0.01 0.9996
 1   rs12497307 MEAN    NA   -0.2 0.54
 1   rs12497307 SD      NA   0.65  1
 1   rs12565286 GENO    C/C   C/G  G/G
 1   rs12565286 COUNTS   3    459  5239
 1   rs12565286 FREQ  0.001  0.08  0.919
 1   rs12565286 MEAN   0.13  -0.02 0.002
 1   rs12565286 SD     0.12  0.96  1.003
 ...

 second input file:
 CHR SNP        VALUE  G11    G12 G22
 1   rs12565286 GENO    C/C   C/G  G/G
 1   rs12565286 COUNTS   1    459  5239
 1   rs12565286 FREQ  0.001  0.28  0.119
 1   rs12565286 MEAN   0.13  -0.05 0.12
 1   rs12565286 SD     0.12  0.9  1.003
 1   rs12497307 GENO   A/A    A/G  G/G
 1   rs12497307 COUNTS  0      2  3000
 1   rs12497307 FREQ    0    0.01 0.996
 1   rs12497307 MEAN    NA   -0.2 0.44
 1   rs12497307 SD      NA   0.65  2
 ...


third input file
 ...


	lambdas <- var.meta(filenames_array=c("cohort1.means", "cohort2.means"))		
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
