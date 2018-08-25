##########################################################################
# Genotyping Uncertainty with Sequencing data - Linkage Disequilibrium (GUSLD)
# Copyright 2017-2018 AgResearch Ltd. <timothy.bilton@agresearch.co.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
#' Deer Dataset
#'
#' Function for extracting the path to the (VCF) file of the deer data set used in the publication
#' by \insertCite{bilton2018genetics2;textual}{GUSbase}.
#'
#' The data consists of 38 SNPs and 704 deer samples genotyped using the genotyping-by-sequencing
#' method. The data is in VCF format (see this \href{https://samtools.github.io/hts-specs/VCFv4.2.pdf}{page}
#' for specification of VCF format).
#'
#' @return A character string of the complete path to the VCF file of the deer data set.
#' @references
#' \insertRef{bilton2018genetics2}{GUSbase}
#' @examples
#' ## extract the name of the vcf file for the deer dataset
#' deerData()
#'
#' @export
## Wrapper function for reading in the Manuka data set for chromosome 11
deerData <- function(){
  return(system.file("extdata", "deer.vcf", package="GUSLD"))
}





