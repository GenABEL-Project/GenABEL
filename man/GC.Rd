\name{GC}
\alias{GC}
\title{Genomic control for various model of inheritance using VIF}
\usage{
  GC(data = 0, p, x, method = "regress", clust = 0,
    vart0 = 0, tmp = 0, min.p = 0.05, CA = FALSE,
    index.filter = 0, n, p.table = 0)
}
\arguments{
  \item{data}{Input vector of Chi square statistic}

  \item{method}{Function of error to be optimized. Can be
  "regress", "median" or "ks.test"}

  \item{p}{Input vector of allele frequencies}

  \item{x}{Model of inheritance (0 for recessive,0.5 for
  additive, 1 for dominant)}

  \item{index.filter}{Indexes for variables that will be
  use for analisis in data vector}

  \item{n}{size of the sample}

  \item{clust}{For developers only}

  \item{vart0}{For developers only}

  \item{tmp}{For developers only}

  \item{min.p}{For developers only}

  \item{CA}{For developers only}

  \item{p.table}{For developers only}
}
\value{
  A list with elements \item{Zx}{output vector corrected
  Chi square statistic} \item{vv}{output vector of VIF}
  \item{exeps}{output vector of exepsons (NA)}
  \item{calrate}{output vector of calrate} \item{F}{F}
  \item{K}{K}
}
\description{
  This function estimates the genomic controls for diffrent
  model (recessive,dominant,additive etc.), using VIF. VIF
  coefficients are estimated by optimizing diffrent error
  functions: regress, median and ks.test.
}
\examples{
#Sorry, this block will appear later
}
\author{
  Yakov Tsepilov
}
\keyword{htest}

