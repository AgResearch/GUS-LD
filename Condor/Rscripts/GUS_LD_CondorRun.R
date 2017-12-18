###############################
# Genotyping Uncertainty with Sequencing data and Linkage Disequilibrium (GUS-LD)
# Copyright 2017 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#########################################################################
#### R script for estimating LD for high-throughput data using approach by Bilton (2017)
#### Script does the calculations over condor and using parallel programming.
#### Script assumes,
####  > Data has been processed in to a genotype matrix (genon) and depth matrix (depth
####  > Correct working directory is specified beforehand
####  > packages 'foreach' and 'doSNOW' can be loaded
####  > Needs to ensure that the file 'GUS_LD_Compute_LDblock.R' is in the right folder.
###
### Author: Timothy P. Bilton
### Date: 16/12/17

### Path to data folder and data file names
dataFolder <- Sys.getenv('DATA_FOLDER')
genonName <- Sys.getenv('GENON_NAME')   # should be a string specifying path
depthName <- Sys.getenv('DEPTH_NAME')   # should be a string specifying path

### Number of blocks the SNPs are divided up into
## Note: The number of queues specified in the condor file needs to be twice this.
nblocks = as.integer(Sys.getenv('NBLOCKS'))

### number of clusters to be used for the parallel program
## Note: need to ensure that the number of cpus asked for in the conder file is equivalent
nClusters <- as.integer(Sys.getenv('NCLUSTERS'))

### Determinine the path from where this file is to the directory where the scripts to compute the LD blocks is located.
path <- Sys.getenv('RSCRIPT_PATH')
runName <- Sys.getenv('RUN_NAME')

### Any additional R libraries required?
RLibrary <- Sys.getenv('R_LIBRARY')


#### Additional, if compute is for specific chromosome
#chrom = gsub(" ", "", chrom, fixed = TRUE)
chrom = ""
if(chrom != ""){ chrom <- paste("_",chrom,sep="")}

### read in the blocks matrix
blocks <- matrix(scan(paste("./",runName,"/LDblocks/BlockMat.txt",sep="")),byrow=T,ncol=4)
block_int <- unique(blocks[,1:2],MARGINS=2)

### Numebr of individuals in the dataset (should be known by the user)
nind <- as.integer(Sys.getenv('NIND'))

### work out which queue we are in 
run <- as.integer(commandArgs(trailingOnly=T))[1] + 1

### select the block we want to work on
g1 <- blocks[run,1:2]
g2 <- blocks[run,3:4]

#### Need to include over R library to access additional packages
.libPaths( c( .libPaths(), RLibrary) )

#### Load the required packages ####
if(!require(foreach) | !require(doSNOW)){
  stop("One of the packages 'foreach' or 'doSNOW' could not be loaded")
}


################### Define additional function required for the computation #################
#### define the likelihood function
pstar.ll <- compiler::cmpfun(function(x,coefMat,id){
  pA1=x[1]; pA2=x[2]; D=x[3]
  pB1 = 1-pA1; pB2 = 1-pA2
  C1 <- max(-pA1*pA2,-pB1*pB2); C2 <- min(pB1*pA2,pA1*pB2)
  if(!((C1 <= D) & (D <= C2))) return(NA)
  else{    
    h1 <- pA1*pA2+D; h2 <- pA1*pB2-D; h3 <- pB1*pA2-D; h4 <- pB1*pB2+D
    p11 <- h1^2
    p21 <- 2*h1*h3
    p31 <- h3^2
    p12 <- 2*h1*h2
    p32 <- 2*h3*h4
    p13 <- h2^2
    p23 <- 2*h2*h4
    p33 <- h4^2
    p22 <- 2*h1*h4+2*h2*h3
    
    return(-(sum(log(p11+(coefMat[,2]*p12+coefMat[,1]*p21+coefMat[,5]*p22)[id[[1]]])) +
               sum(log((coefMat[,3]*p21+coefMat[,3]*coefMat[,2]*p22)[id[[2]]])) +
               sum(log(p31 + (coefMat[,1]*p21+coefMat[,2]*p32+coefMat[,5]*p22)[id[[3]]])) +
               sum(log((coefMat[,4]*p12+coefMat[,4]*coefMat[,1]*p22)[id[[4]]])) +
               sum(log((coefMat[,3]*coefMat[,4]*p22)[id[[5]]])) + 
               sum(log((coefMat[,4]*p32+coefMat[,4]*coefMat[,1]*p22)[id[[6]]])) +
               sum(log(p13 + (coefMat[,2]*p12+coefMat[,1]*p23+coefMat[,5]*p22)[id[[7]]])) +
               sum(log((coefMat[,3]*p23+coefMat[,3]*coefMat[,2]*p22)[id[[8]]])) +
               sum(log(p33 + (coefMat[,2]*p32+coefMat[,1]*p23+coefMat[,5]*p22)[id[[9]]]))
    ))
  }
})

### function for estimating LD estimates for given MLE values
LDest <- function(MLE,n,measure){
  D <- MLE[3]#*(2*n)/(2*n-1)
  C1 = max(-prod(MLE[1:2]),-prod(1-MLE[1:2])); C2 = min((1-MLE[1])*MLE[2],MLE[1]*(1-MLE[2]))
  D <- ifelse(D<0,max(D,C1),min(D,C2))
  if(measure=='Dprime'){
    Dprime <- D/ifelse(D<0,C1,C2)
    return(Dprime)
  } else if(measure=='r2'){
    r2 <- D^2/(prod(c(MLE[1:2],1-MLE[1:2])))
    return(r2)
  }
}

### Function to define how to put the vecotr back together for the foreach loop
comb <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

####################### Now look to do the computations ####################
#### set up the clusters
attempts <- 0
while(attempts<10){
  attempts= attempts+1
  madeClus <- try(cl <- makeCluster(nClusters)) # define the number of clusters
  if(all(class(madeClus)=="try-error")) next
  madeClus <- try(registerDoSNOW(cl)) # required or else foreach will just run sequentially 
  if(class(madeClus)!="try-error") break
  else stopCluster(cl)
}

computedLD <- try(source(paste(path,'/GUS_LD_ComputeLDblock.R',sep=""))) # try command will ensure this script will keep going even if there is an error with the sourced script

if(exists(deparse(substitute(cl))))
  stopCluster(cl) # essential to do this if we set up the clusters

## Check to see if the LD block was actually computed, if not throw an error and stop the script
if(class(computedLD)=='try-error'){
  stop("Problem with computing the LD block matrix. Check the error file (.err) in the log folder to find the cause.")
}

#### Now want to write the matrix to a file.
## Note that we are putting them into the a subfolder.
write.table(newblock[[1]],file=paste("./",runName,"/LDblocks/Dprime_Block",run,chrom,".txt",sep=""))
write.table(newblock[[2]],file=paste("./",runName,"/LDblocks/r2_Block",run,chrom,".txt",sep=""))








