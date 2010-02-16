\name{phdata}
\title{phdata}
\usage{phdata(data)
gtdata(data) 
chromosome(data)
coding(data)
idnames(data)
map(data)
nids(data)
nsnps(data)
snpnames(data)
strand(data)}
\description{Access to slots of \link{gwaa.data-class} object}
\details{Functions to access slots of \link{gwaa.data-class} object. In particular:

phdata: phenotypic data

gtdata: genotypic data

idnames: subject IDs

snpnames: marker IDs

nids: number of subjects

nsnps: number of markers

...}
\alias{gtdata}
\alias{chromosome}
\alias{coding}
\alias{gtdata}
\alias{idnames}
\alias{map}
\alias{nids}
\alias{nsnps}
\alias{phdata}
\alias{snpnames}
\alias{strand}
\value{content of the particular slot}
\arguments{\item{data}{object of \link{gwaa.data-class}}}
