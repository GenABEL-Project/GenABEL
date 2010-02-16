\name{mach2databel}
\alias{mach2databel}
\title{mach2databel}
\usage{mach2databel(mldosefile, mlinfofile, outfile)}
\description{converts MACH-imputed files to DatABEL (filevector) format}
\details{This function converts mach-imputed files to \code{DatABEL} (filevector) format.
After conversion, two files (outfile.fvi and outfile.fvd), corresponding 
to single filevector object, will appear on the disk; databel_filtered_R 
object connected to these files will be returned to R}
\value{databel_filtered_R-class object}
\author{Yurii Aulchenko}
\keyword{IO}
\keyword{manip}
\arguments{\item{mldosefile}{MACH mldose file name}
\item{mlinfofile}{MACH mlinfo file name}
\item{outfile}{output file name}}
