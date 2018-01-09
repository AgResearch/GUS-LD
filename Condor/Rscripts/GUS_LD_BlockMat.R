###############################
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
#### File putting putting groups of SNP's into blocks
#### based on the specifed number of blocks of the data matrix.
#### It then creates the required files to compute all the LD pairwise 
#### calculations using condor
####
###  Author: Timothy P. Bilton
###  Date: 16/12/17

### extract the environment variables needed.
### Number of blocks the SNPs are divided up into
nblocks = as.integer(Sys.getenv('NBLOCKS'))

### number of clusters to be used for the parallel program
nClusters <- as.integer(Sys.getenv('NCLUSTERS'))

### Path to data folder and data file names
dataFolder <- Sys.getenv('DATA_FOLDER')
genonName <- Sys.getenv('GENON_NAME')   # should be a string specifying path
depthName <- Sys.getenv('DEPTH_NAME')   # should be a string specifying path

### Name of this particular run.
runName <- Sys.getenv('RUN_NAME')

### Any additional R libraries required?
RLibrary <- Sys.getenv('R_LIBRARY')

#### Load data
genon <- as.matrix(read.table(paste(dataFolder,"/",genonName,sep=""),comment.char="",colClasses=integer(0)))
depth <- as.matrix(read.table(paste(dataFolder,"/",depthName,sep=""),comment.char="",colClasses=integer(0)))

## need to make sure that the maximum depth is 1000 (for numeric reasons so that the optimization works)
depth[which(depth > 500)] <- 500

totalSnps <- ncol(genon) # Total number of SNPs in the matrix
nind <- nrow(genon)      # Total number of individuals

### determine the blocks of the matrix
cuts <- floor(seq(0,totalSnps,length.out=nblocks+1))
block_int <- cbind(cuts[-(nblocks+1)]+1,cuts[-1])
blocks <- matrix(c(block_int[rep(1:nblocks,nblocks:1),],
                   block_int[sequence(nblocks:1) + rep(0:nblocks,nblocks:0),]),ncol=4)

## write this block matrix to a file (since we will need it later on)
write.table(blocks,paste("./",runName,"/LDblocks/BlockMat.txt",sep=""),row.names=F,col.names=F)

### Save the data as spearate text files in the temporary data directory
for(block in 1:nrow(block_int)){
  write.table(genon[,block_int[block,1]:block_int[block,2]],file=paste("./",runName,"/tempData/genon",block,".txt",sep=""),
              row.names=F,col.names=F)
  write.table(depth[,block_int[block,1]:block_int[block,2]],file=paste("./",runName,"/tempData/depth",block,".txt",sep=""),
              row.names=F,col.names=F)
}

## additional environment variable needed for condor file
path <- Sys.getenv('RSCRIPT_PATH')

### work out a reasonable memory requrest size.
memory <-  Sys.getenv('CONDOR_MEMORY')

## Now generate the condor file...

cat('Executable     = ',runName,'/',runName,'.sh
Universe        = Vanilla
arguments       = $(Process)
environment     = NBLOCKS=',nblocks,';NCLUSTERS=',nClusters,';GENON_NAME=',genonName,';DEPTH_NAME=',depthName,';NIND=',nind,';RSCRIPT_PATH=',path,';DATA_FOLDER=',dataFolder,';RUN_NAME=',runName,';R_LIBRARY=',RLibrary,'
request_cpus    = ',nClusters,'
request_memory  = ',memory,'
error           = ./',runName,'/log/condor_',runName,'_run$(Process).err
output          = ./',runName,'/log/condor_',runName,'_run$(Process).out
log             = ./',runName,'/log/condor_',runName,'_run$(Process).log
queue ',sum(seq(1:nblocks)), file=paste(runName,"/",runName,".condor",sep=""),sep="")

### write the condor Executable file.
cat('#!/bin/bash

/usr/bin/Rscript ',path,'/GUS_LD_CondorRun.R $1', file=paste(runName,"/",runName,".sh",sep=""),sep="")

