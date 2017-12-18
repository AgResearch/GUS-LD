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

echo "<<<<< remove environment variables created >>>>>"
unset NBLOCKS
unset NCLUSTERS
unset GENON_NAME
unset DEPTH_NAME
unset DATA_FOLDER
unset CONDOR_MEMORY

echo "<<<< Construct the full LD matrix >>>>"
/usr/bin/Rscript $RSCRIPT_PATH$'/GUS_LD_CombineBlocks.R' $1

# The resulting matrix of LD values will be saved as the 
# RData file. To read into R use the command
echo "< To load the data into R, use the code >"
echo
echo "load('"$RUN_NAME"/"$RUN_NAME"_LDmat.RData')"
echo
echo "< provided that the current directory in R is LDmapping >" 
echo
echo "< Creating Heatmaps>"

# Tidy-up, remove the data created so that we are not taking up too much disk space
# rm -rf $RUN_NAME/tempData

### remove the remaining environment variables
unset RSCRIPT_PATH
unset RUN_NAME
