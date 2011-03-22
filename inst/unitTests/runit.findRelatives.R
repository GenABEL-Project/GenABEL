### --- Test setup ---

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library(RUnit)
    library(GenABEL)
    library(DatABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.findRelatives <- function()
{
    data(ge03d2.clean)
    nloci <- 2000
    set.seed(1)
    df <- ge03d2.clean[,sort(sample(autosomal(ge03d2.clean),nloci))]
    q <- summary(gtdata(df))$"Q.2"
#
# g1---g2
#    |
#    g3----g4
#       /\
#     g5  g6---g7
#            |
#            g8---g9
#               |
#               g10
    nids <- 10
    founders <- sample(idnames(df),5)
    gt <- matrix(ncol=nloci,nrow=nids)
    gt[1,] <- rbinom(nloci,2,q)
    gt[2,] <- rbinom(nloci,2,q)
    gt[4,] <- rbinom(nloci,2,q)
    gt[7,] <- rbinom(nloci,2,q)
    gt[9,] <- rbinom(nloci,2,q)
    gt[3,] <- generateOffspring(gt[1,],gt[2,],q=q)
    gt[5,] <- generateOffspring(gt[3,],gt[4,],q=q)
    gt[6,] <- generateOffspring(gt[3,],gt[4,],q=q)
    gt[8,] <- generateOffspring(gt[6,],gt[7,],q=q)
    gt[10,] <- generateOffspring(gt[8,],gt[9,],q=q)
    a<-findRelatives(gt,q=q,nmei=c(1:2))
	checkIdentical(a$compressedGuess[1,3],"1")
	checkIdentical(a$compressedGuess[2,3],"1")
	checkIdentical(a$compressedGuess[1,5],"2")
	checkIdentical(a$compressedGuess[1,6],"2")
	checkIdentical(a$compressedGuess[2,5],"2")
	checkIdentical(a$compressedGuess[2,6],"2")
	checkIdentical(a$compressedGuess[5,6],"2+2")
}