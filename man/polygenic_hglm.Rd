\name{polygenic_hglm}
\alias{polygenic_hglm}
\title{Estimation of polygenic model...}
\usage{polygenic_hglm(formula, kinship.matrix, data, family=gaussian(),
    conv=1e-08, maxit=100, ...)}
\description{Estimation of polygenic model}
\details{Estimation of polygenic model using hierarchical
GLM (hglm package)
Tell when REML is used, what otherwise; 
what it is doing, how to do tests, add more examples, references, etc.}
\author{Xia Shen, Yurii Aulchenko}
\references{need to put reference here}
\seealso{\code{\link{polygenic}},
\code{\link{mmscore}},
\code{\link{grammar}}}
\keyword{htest}
\arguments{\item{formula}{Formula describing fixed effects to be used in analysis, e.g. 
y ~ a + b means that outcome (y) depends on two covariates, a and b. 
If no covariates used in analysis, skip the right-hand side of the 
equation.}
\item{kinship.matrix}{Kinship matrix, as provided by e.g. ibs(,weight="freq"), 
or estimated outside of GenABEL from pedigree data.}
\item{data}{An (optional) object of \code{\link{gwaa.data-class}} or a data frame with 
outcome and covariates}
\item{family}{a description of the error distribution and link function to be 
used in the mean part of the model (see \code{\link{family}} for details of 
family functions)}
\item{conv}{'conv' parameter of 'hglm' (stricter then default)}
\item{maxit}{'maxit' parameter of 'hglm' (stricter then default)}
\item{...}{other parameters passed to 'hglm' call}}
\examples{data(ge03d2ex.clean)
df <- ge03d2ex.clean[,autosomal(ge03d2ex.clean)]
gkin <- ibs(df,w="freq")

# for quantitative traits
h2ht <- polygenic_hglm(height ~ sex + age,kin=gkin,df)
# estimate of heritability
h2ht$esth2
# other parameters
h2ht$h2an

# for binary traits
h2dm <- polygenic_hglm(dm2 ~ sex + age,kin=gkin,df,family=binomial(link = 'logit'))
# estimated parameters
h2dm$h2an}
