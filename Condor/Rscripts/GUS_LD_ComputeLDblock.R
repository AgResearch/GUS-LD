#########################################################################
# Genotyping Uncertainty with Sequencing data and Linkage Disequilibrium (GUS-LD)
# Copyright (C) 2017 AgResearch Ltd.
#
# GUS-LD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GUS-LD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#########################################################################
#### File for computing all pairwise LD for a given block, whether diagonal or off diagonal,
#### using parallel programming in R.
### Author Timothy P. Bilton
### Date: 16/12/17


### Define inportant variables
my.grid <- lapply(1:9, function(i) as.matrix(expand.grid(2:0,2:0))[i,])
nSnps <- g1[2]-(g1[1]-1)
Dprime <- as.integer(numeric(nSnps))
r2 <- as.integer(numeric(nSnps))

#### For the diagonal block
if(all(g1==g2)){
  # Load the data
  blockNo <- which(apply(block_int,1,function(x) all(x==g1)))
  newGenon <- matrix(scan(paste(runName,"/tempData/genon",blockNo,".txt",sep=""),
                            what=integer(0)),nrow=nind,ncol=g1[2]-(g1[1]-1),byrow=T)
  newDepth <- matrix(scan(paste(runName,"/tempData/depth",blockNo,".txt",sep=""),
                          what=integer(0)),nrow=nind,ncol=g1[2]-(g1[1]-1),byrow=T)
  
  # run the foreach loop to estimate all the pairwise LD values
  newblock <- foreach(rowSnp=iter(1:nSnps), .combine='comb',.multicombine=TRUE) %dopar% {
    for(colSnp in seq_len(rowSnp-1)){
      ind <- c(rowSnp,colSnp)
      GBSgeno <- newGenon[,ind]
      obs <- complete.cases(GBSgeno)
      initialV <- c(0.5,0.5,0)
      # outline the variables
      GBSgeno <- matrix(GBSgeno[obs,],ncol=2)
      depth2 <- matrix(newDepth[obs,ind], ncol=2)
      
      K <- 1/(2^depth2)
      coefMat <- cbind(K,1-2*K,K[,1]*K[,2])
      n <- sum(obs) # number of observations
      id <- lapply(my.grid,function(x) which(GBSgeno[,1]==x[1] & GBSgeno[,2]==x[2]))
      
      ## compute MLE of pA1, pA2 and D
      MLE <- try(optim(initialV,pstar.ll,method="Nelder-Mead", control=list(reltol=1e-15),
                   coefMat=coefMat,id=id)$par)
      if(class(MLE)=='try-error'){
        Dprime[colSnp] <- NA
        r2[colSnp] <- NA
      }else{
        Dprime[colSnp] <- LDest(MLE,n,measure='Dprime')
        r2[colSnp] <- LDest(MLE,n,measure='r2')
      }
    }
    return(list(as.integer(round(Dprime*100)),as.integer(round(r2*100))))
  }
} else{ # for the non-diagonal blocks
  ## load data into 2 parts
  # part 1: row snps
  blockNo1 <- which(apply(block_int,1,function(x) all(x==g1)))
  newGenon <- matrix(scan(paste(runName,"/tempData/genon",blockNo1,".txt",sep=""),
                          what=integer(0)),nrow=nind,ncol=g1[2]-(g1[1]-1),byrow=T)
  newDepth <- matrix(scan(paste(runName,"/tempData/depth",blockNo1,".txt",sep=""),
                          what=integer(0)),nrow=nind,ncol=g1[2]-(g1[1]-1), byrow=T)
  
  # part 2: col snps
  blockNo2 <- which(apply(block_int,1,function(x) all(x==g2)))
  newGenon2 <- matrix(scan(paste(runName,"/tempData/genon",blockNo2,".txt",sep=""),
                          what=integer(0)),nrow=nind,ncol=g2[2]-(g2[1]-1),byrow=T)
  newDepth2 <- matrix(scan(paste(runName,"/tempData/depth",blockNo2,".txt",sep=""),
                          what=integer(0)),nrow=nind,ncol=g2[2]-(g2[1]-1), byrow=T)
  
  rowLen <- g2[2]-(g2[1]-1)
  # run the foreach loop to estimate all the pairwise LD values
  newblock <- foreach(rowSnp=iter(1:rowLen), .combine='comb',.multicombine=TRUE) %dopar% {
    extraGenon <- newGenon2[,rowSnp]; extraDepth <- newDepth2[,rowSnp]
    freq <- c(colMeans(newGenon,na.rm=T)/2,mean(extraGenon,na.rm=T)/2)
    for(colSnp in 1:nSnps){
      GBSgeno <- cbind(newGenon[,colSnp],extraGenon)
      obs <- complete.cases(GBSgeno)
      initialV <- c(0.5,0.5,0) #c(freq[c(colSnp,nSnps+1)],0)
      # outline the variables
      GBSgeno <- matrix(GBSgeno[obs,],ncol=2)
      depth2 <- cbind(newDepth[obs,colSnp],extraDepth[obs])
            
      K <- 1/(2^depth2)
      coefMat <- cbind(K,1-2*K,K[,1]*K[,2])
      n <- sum(obs) # number of observations
      id <- lapply(my.grid,function(x) which(GBSgeno[,1]==x[1] & GBSgeno[,2]==x[2]))
      
      ## compute MLE of pA1, pA2 and D
      MLE <- try(optim(initialV,pstar.ll,method="Nelder-Mead", control=list(reltol=1e-15),
                   coefMat=coefMat,id=id)$par)
      if(class(MLE)=='try-error'){
        Dprime[colSnp] <- NA
        r2[colSnp] <- NA
      }else{
        Dprime[colSnp] <- LDest(MLE,n,measure='Dprime')
        r2[colSnp] <- LDest(MLE,n,measure='r2')
      }
    }
    return(list(as.integer(round(Dprime*100)),as.integer(round(r2*100))))
  }
}

