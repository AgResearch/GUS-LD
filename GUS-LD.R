#########################################################################
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
##### File for computing LD estimates using high-throughtput data
### Author: Timothy P. Bilton
### Date: 16/12/17

######### Notes ###########
# We require:
#   > the data set is coded in terms of 0,1 and 2
#     representing the number of alleles present for the major allele
#   > The depths are gives 


###### Define functions

## likelihood function
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

## function for estimating LD estimates for given MLE values
Iden_D <- function(pA1,pA2,D){
  return(D)
}
Dprime <- function(pA1,pA2,D){
  C1 = max(-prod(c(pA1,pA2)),-prod(1-c(pA1,pA2))); C2 = min((1-pA1)*pA2,pA1*(1-pA2))
  return(D/ifelse(D<0,C1,C2))
}
r2 <- function(pA1,pA2,D){
  return(D^2/(prod(c(pA1,pA2,1-c(pA1,pA2)))))
}

## needed for foreach loop
comb <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

############ Main Function ########################
#### Computing matrics of the GBS Dprime and r2 estimates
# Arguments are
#  - GBSdata: Matrix of GBS genotype calls. 0=BB, 1=AB, 2=AA, NA=missing
#  - depths: Matrix of read depths associated with the GBS genotype calls
#  - parallel: If TRUE then parallel computing is used to speed up the computations
#  - Nclust: Number of clusters the foreach loop using for the parallel computing
#            Note: make sure you don't use for than the number of cores available
#  - LDmeasure: A vector or function numbers which compute the required LD measure from the MLE of D
#               Note: Some predefined functions are above.
GUS_LD <- function(genon,depths,Nclust=4,parallel=TRUE,LDmeasure=c('Dprime','r2')){
  
  ## Do some checks
  if(!is.matrix(genon)|!is.matrix(depths))
    stop("Matrix objects required for inputs")
  if(any(!(na.omit(unique(as.vector(genon))) %in% c(0,1,2))))
    stop("GBS data has invalid input. Only 0, 1, 2 and NA's are allowed")
  if(any(na.omit(as.vector(depths)<0)))
    stop("Negative depths are invalid.")
  if(any(dim(genon)!=dim(depths)))
    stop("GBS data and depths matrices are different dimensions")
  if(any(!is.na(genon[which(depths==0)])))
    stop("Non missing GBS calls with zero depth")
  if(!(length(LDmeasure) > 0))
    stop("No LD measure has been specified")
  if(Nclust < 0){
    warning("Number of clusters specified in the parallelization is not positive.\nUsing no paralelization instead.")
    parallel = FALSE  
  }
  # If to run it using parallel programming
  if(parallel){
    #### Load the required packages ####
    if(!require(foreach) | !require(doSNOW)){
      stop("One of the packages 'foreach' or 'doSNOW' could not be loaded")
    }
    
    cl <- makeCluster(Nclust) # define the number of clusters
    registerDoSNOW(cl) # required or else foreach will just run sequentially 
    
    freq <- colMeans(genon,na.rm=T)/2
    nSnps <- ncol(genon)
    r2.ms <- matrix(nrow=nSnps,ncol=nSnps)
    my.grid <- lapply(1:9, function(i) as.matrix(expand.grid(2:0,2:0))[i,])
    
    status <- try(LDmat <- foreach(snp1=iter(1:nSnps),.combine=comb,.multicombine=TRUE,
                                   .export=c("pstar.ll",LDmeasure)) %dopar% {
                                     LDvec <- vector(mode='list')
                                     LDvec <- replicate(length(LDmeasure),numeric(nSnps),simplify=F)
                                     for(snp2 in seq_len(snp1-1)){
                                       ind <- c(snp1,snp2)
                                       GBSgeno <- cbind(genon[,ind])
                                       obs <- complete.cases(GBSgeno)
                                       initialV <- c(freq[ind],0)
                                       
                                       # outline the variables
                                       GBSgeno <- matrix(GBSgeno[obs,],ncol=2)
                                       depth2 <- matrix(depths[obs,ind], ncol=2,byrow=F)
                                       
                                       K <- 1/(2^depth2)
                                       coefMat <- cbind(K,1-2*K,K[,1]*K[,2])
                                       n <- sum(obs) # number of observations
                                       id <- lapply(my.grid,function(x) which(GBSgeno[,1]==x[1] & GBSgeno[,2]==x[2]))
                                       
                                       # Determine if homozygous is with respect to the major or minor allele
                                       ## compute MLE of pA1, pA2 and D
                                       MLE <- optim(par=initialV,fn=pstar.ll,method="Nelder-Mead", control=list(reltol=1e-15),coefMat=coefMat,id=id)$par
                                       D <- MLE[3]*(2*n)/(2*n-1)
                                       C1 = max(-prod(MLE[1:2]),-prod(1-MLE[1:2])); C2 = min((1-MLE[1])*MLE[2],MLE[1]*(1-MLE[2]))
                                       D <- ifelse(D<0,max(D,C1),min(D,C2))
                                       for(meas in 1:length(LDmeasure)) {
                                         LDvec[[meas]][snp2] <- get(LDmeasure[[meas]])(pA1=MLE[1],pA2=MLE[2],D=D)
                                       }
                                     }
                                     return(lapply(LDvec,round,digits=4))
                                   } )
    
    if(class(status)=='try-error'){
      stopCluster(cl)
      stop('Error: Problem with estimation proceedure using parallel computing.')
    } else {  stopCluster(cl) }
    
    names(LDmat) <- LDmeasure
    
    for(i in 1:length(LDmeasure)){
      diag(LDmat[[i]]) <- 1
      LDmat[[i]][upper.tri(LDmat[[i]])] <- t(LDmat[[i]])[upper.tri(t(LDmat[[i]]))]
    }
    
    return(LDmat)
  }
  else{
    freq <- colMeans(genon,na.rm=T)/2
    nSnps <- ncol(genon)
    LDmat <- replicate(length(LDmeasure),matrix(nrow=nSnps,ncol=nSnps),simplify=FALSE)
    names(LDmat) <- LDmeasure
    my.grid <- lapply(1:9, function(i) as.matrix(expand.grid(2:0,2:0))[i,])
    
    for(i in 1:(nSnps-1)){
      for(j in (i+1):nSnps){
        ind <- c(i,j)
        obs <- complete.cases(genon[,ind])
        initialV <- c(freq[ind],0)
        # outline the variables
        GBSgeno <- matrix(genon[obs,ind],ncol=2)
        depth <- depths[obs,ind]
        
        K <- matrix(1/(2^depth), ncol=2)
        coefMat <- cbind(K,1-2*K,K[,1]*K[,2])
        n <- sum(obs) # number of observations
        id <- lapply(my.grid,function(x) which(GBSgeno[,1]==x[1] & GBSgeno[,2]==x[2]))
        
        ## compute MLE of pA1, pA2 and D
        MLE <- optim(par=initialV,fn=pstar.ll,method="Nelder-Mead", control=list(reltol=1e-15),coefMat=coefMat,id=id)$par
        D <- MLE[3]*(2*n)/(2*n-1)
        C1 = max(-prod(MLE[1:2]),-prod(1-MLE[1:2])); C2 = min((1-MLE[1])*MLE[2],MLE[1]*(1-MLE[2]))
        D <- ifelse(D<0,max(D,C1),min(D,C2))
        for(meas in 1:length(LDmeasure)) {
          LDmat[[meas]][i,j] <- get(LDmeasure[[meas]])(pA1=MLE[1],pA2=MLE[2],D=D)
        }
      }
    }
    for(i in 1:length(LDmeasure)){
      diag(LDmat[[i]]) <- 1
      LDmat[[i]][lower.tri(LDmat[[i]])] <- t(LDmat[[i]])[lower.tri(t(LDmat[[i]]))]
    }
    return(LDmat)    
  }
}




