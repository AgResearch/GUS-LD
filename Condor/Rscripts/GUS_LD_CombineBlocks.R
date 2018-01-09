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
#### R script that takes the LD blocks generated using the 'GUS_LD_CondorRun.R',
#### combines them into a half matrix and saves to an RData file.
####
### Author:  Timothy P. Bilton  
### Date: 16/12/17

## Access to any additional R library
RLibrary <- Sys.getenv('R_LIBRARY')
.libPaths( c( .libPaths(), RLibrary) )

runName <- Sys.getenv('RUN_NAME')

### Function that takes a list of LDblocks and combines them into a matrix. 
### Note: It assumes that the blocks are to be added sequentially down a column
combineBlocks <- function(List,blocks){
  nb <- nrow(unique(blocks[,1:2]))
  bLen <- blocks[1:nb,4]-blocks[1:nb,3]+1
  
  LDmat <- do.call(rbind,List[c(1:nb)])
  for(i in seq(from=1,length.out=(nb-1))){
    addBlock <- matrix(0,nrow=sum(bLen[1:i]),ncol=bLen[i+1])
    LDmat <- cbind(LDmat,do.call(rbind,c(list(addBlock),List[seq(from=sum(nb+1-1:i)+1,length.out=nb-i)])))
  }
  return(LDmat)
}


#### First generate a list that contains all the subblocks that we need and combine to get a matrix.
for(LDmeasure in c('Dprime','r2')){

  blockList <- list() #list to contain all the LDblocks
  ## read in the block matrix
  blocks <- matrix(scan(paste("./",runName,"/LDblocks/BlockMat.txt",sep="")),byrow=T,ncol=4)
  
  # Now out the blocks into a list.
  for(blockNo in 1:nrow(blocks)){
    blockList <- c(blockList,list(as.matrix(read.table(paste("./",runName,"/LDblocks/",LDmeasure,"_Block",blockNo,".txt",sep="")))))
  }
  
  # combine the blocks into a matrix
  LDmat <- combineBlocks(blockList,blocks)
  
  ### If required. Turn the half matrix into a full symmetric matrix 
  diag(LDmat) <- 100
  LDmat[upper.tri(LDmat)] <- t(LDmat)[upper.tri(t(LDmat))]
  
  ### covert data to actural floats instead of integers and save an RData file
  assign(LDmeasure,unname(LDmat)/100)
}

## save to an R workspace
save('Dprime','r2',file=paste(runName,"/",runName,"_LDmat.RData",sep=""))





