## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set( collapse = TRUE, comment = "#>", eval=TRUE, class.source="courier-font")
library(GUSLD)

## ----simDS---------------------------------------------------------------
vcffile <- deerData() # extract filename
vcffile               # filename stored in object vcffile

## ----loadData------------------------------------------------------------
# convert VCF file to an RA file
rafile <- VCFtoRA(infilename = vcffile, direct = "./")

## ----rafile--------------------------------------------------------------
rafile # file path of RA file

## ----readRA--------------------------------------------------------------
RAdata <- readRA(rafile = rafile, sampthres = 0.01, excsamp = NULL)

## ----class_RAdata--------------------------------------------------------
class(RAdata)

## ----RAdata, class.source="courier-font"---------------------------------
RAdata

## ----makeUR--------------------------------------------------------------
urpop <- makeUR(RAobj = RAdata, filter = list(MAF = 0.01, MISS = 0.5, HW = c(-0.05, Inf), 
                                              MAXDEPTH = 500), nThreads = 2)

## ----urpop---------------------------------------------------------------
urpop

## ----allSNPs-------------------------------------------------------------
LDres <- GUSLD(URobj = urpop, nClust = 2)

## ----allSNPs_output------------------------------------------------------
head(LDres)

## ----allSNPs_file--------------------------------------------------------
LDres <- GUSLD(urpop, filename="LD_deer", dp = 4)

## ----subet---------------------------------------------------------------
## set up the SNP pairs
pairs <- as.matrix(expand.grid(1:10,11:20))
head(pairs)

## ----LDres_block---------------------------------------------------------
## Compute the LD estimates between SNPs in differen blocks
LDres <- GUSLD(urpop, SNPpairs = pairs)
head(LDres)

## ----compute_LD----------------------------------------------------------
# To view the function for D'
Dprime
# To compute D' for D=0.1, pA1=0.5, pA2=0.75
Dprime(pA1=0.5, pA2=0.75, D=0.1)

## ----multiple_LDmeasure--------------------------------------------------
LDres <- GUSLD(urpop, LDmeasure=c("Dcoef","Dprime","r2"))
head(LDres)

## ----LD_r----------------------------------------------------------------
## Define the correlation coefficient
LD_r <- function(pA1,pA2,D){
return( D/sqrt((prod(c(pA1,pA2,1-c(pA1,pA2))))) )
}

## ----LDmes_r-------------------------------------------------------------
LDres <- GUSLD(urpop, LDmeasure = "LD_r")
head(LDres)

