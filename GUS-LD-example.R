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
#### Script giving example of how to use GUS-LD 
### Author: Timothy P. Bilton
### Date: 16/12/17

## Read in the data
load("deer.RData")

## Source the function for calculating the pairwise LD
source("GUS-LD.R")

## Compute the LD matrix
LD_mat <- GUS_LD(genon=genon, depth_Ref = depth_Ref, depth_Alt = depth_Alt, parallel = F)

## Plot the results
# D'
image(LD_mat[[1]],col=rev(heat.colors(100)),main="Plot of D'")
# r2
image(LD_mat[[2]],col=rev(heat.colors(100)),main=expression("Plot of"~r^2))

### Can use parallelization to speed up computation.
### Note: Nclust gives tehe number of clusters used in the parallelization
###       DO NOT set as being equal to or larger than the number of cores in your computer. 
LD_mat <- GUS_LD(genon=genon, depth_Ref = depth_Ref, depth_Alt = depth_Alt, Nclust = 3, parallel = T)

### Can specify an alternative LD measure
##  Measures already defined
##   - Iden_D: D  = pA1pA2 - D
##   - Dprime: D' = D/max(D)
##   - r2: r2 = D^2/(pA1(1-pA1)pA2(1-pA2))

## Define the correlation coefficient LD measure
LD_r <- function(pA1,pA2,D){
  return( D/sqrt((prod(c(pA1,pA2,1-c(pA1,pA2))))) )
}

# Compute the correlation coefficient
LD_mat2 <- GUS_LD(genon=genon, depth_Ref = depth_Ref, depth_Alt = depth_Alt, Nclust = 3, parallel = T, LDmeasure = "LD_r")

image(LD_mat2$LD_r,col=rev(heat.colors(100)),main=expression("Plot of r"))

## To extract the error parameters
image(LD_mat2$epsilon,col=rev(heat.colors(100)),main=expression("Plot of sequencing error"))
hist(LD_mat2$epsilon)

#### Code for reading in an VCFfile
source("readVCF.R")

### Arguments for readVCF function:
##  - genofile: Nmae of VCFfile
##  - gform: If aligned to reference, use reference. If called using UNEAK, then use 'uneak'.
##  - sampthres: Minimum sample depth at which samples are to be removed.
library(RCurl)
filename = "https://raw.github.com/tpbilton/GUSMap/master/inst/extdata/Manuka_chr11.vcf"
newData <- readVCF(infilename=getURL(filename), gform="reference")

LD_mat2 <- GUS_LD(genon=newData$genon[,1:30], depth_Ref = newData$depth_Ref[,1:30], depth_Alt = newData$depth_Alt[,1:30],
                  Nclust = 3, parallel = T)

image(LD_mat2$r2,col=rev(heat.colors(100)),main=expression("Plot of r"))

## if SNPs where called using uneak. Then set gform="uneak"




