\name{PGC}
\alias{PGC}
\title{Polynomial genomic control}
\usage{
  PGC(data, method = "regress", p, df, pol.d = 3,
    plot = TRUE, index.filter = 0, start.corr = FALSE)
}
\arguments{
  \item{data}{Input vector of Chi square statistic}

  \item{method}{Function of error to be optimized. Can be
  "regress", "median" or "ks.test"}

  \item{p}{Input vector of allele frequencies}

  \item{df}{Number of dergees of freedom}

  \item{pol.d}{Polinom degree}

  \item{plot}{If true, functon makes plot of lambda from
  allele frequencies}

  \item{start.corr}{For regress method use it only when you
  want to make calculations faster}

  \item{index.filter}{Index of variables in data vector,
  that will be used in analysis, if zero - all variables
  will be used}
}
\value{
  A list with values \item {data} output vector corrected
  Chi square statistic \item {b} polinom coefficients
}
\description{
  This function estimates the genomic controls for diffrent
  model and degrees of freedom, using polinom function.
  Polinomic coefficients are estimated by optimizing
  diffrent error functions: regress, median and ks.test.
}
\examples{
a=1
b=1
a+b
}
\author{
  Yakov Tsepilov
}
\keyword{'}
\keyword{htest}

