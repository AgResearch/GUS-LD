#!/bin/bash

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

echo "<<<<< Create environment variables required to run program >>>>>"
# Define path to data folder and data file names
export GENON_NAME=
export DEPTH_NAME=
export DATA_FOLDER=

# specify the name of this run
export RUN_NAME=

# Define the number of blocks to split the LD matrix into
export NBLOCKS=20

# Define the number of clusters to use in foreach loop
# Four has been a good number
export NCLUSTERS=4

# specify the amount of memory required for each condor job.
# Needs to be montored!!!!
# Need more for larger block sizes and for more clusters in the parallel processing
export CONDOR_MEMORY=500

# specify any additional paths to R libraries that my be required
export R_LIBRARY=

echo "<<<<< Creating additional environment variables and folders >>>>>"
# Are not to be changed (provided this bash file is called from LDmapping folder 
# specify path to R scripts
export RSCRIPT_PATH="./Rscripts"

# make the required directories 
mkdir $RUN_NAME  
mkdir $RUN_NAME/LDblocks
mkdir $RUN_NAME/log 

# make folder for the tempory datasets
mkdir $RUN_NAME/tempData


echo "<<<< Generating the blocks for the LD matrix >>>>"
/usr/bin/Rscript $RSCRIPT_PATH$'/GUS_LD_BlockMat.R' $1

# convert the created sh file to linux format
dos2unix $RUN_NAME/$RUN_NAME.sh

# need to allow access to the above file
chmod +x $RUN_NAME/$RUN_NAME.sh

echo "<<<< Sumitting jobs to condor >>>>"
# once the above R script has run, want to submit the jobs on condor
condor_submit $RUN_NAME/$RUN_NAME.condor

echo "<<<< condor jobs submitted >>>>"

## Run the script GUS_LD_ComputeMat.sh after the conder has finished running to produce the LD matrix.



