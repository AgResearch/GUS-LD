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
#' Linkage Disequilibrium (LD) Measures
#'
#' Functions for computing LD measures for predefined (or estimated) allele frequency values and
#' a predefined (or estimated) disequilibrium coefficient value.
#'
#' The LD measures available in GUSLD are:
#' \describe{
#' \item{\code{Dcoef}}{The disequilibrium coefficient \insertCite{lewontin1960evolution}{GUSLD}
#' defined as \eqn{D=p_{A_1A_2}-p_{A_1}p_{A_2}}}
#' \item{\code{Dprime}}{The normalized disequilibrium coefficient \insertCite{lewontin1964genetics}{GUSLD} defined as
#' \eqn{D'=D/D_{max}} where \eqn{D_{max} = min(p_{A_1}p_{A_2},(1-p_{A_1})(1-p_{A_2}))} if \eqn{D>0} and
#' \eqn{D_{max} = min(p_{A_1}(1-p_{A_2}),(1-p_{A_1})p_{A_2})}}
#' \item{\code{r2}}{The squared correlation coefficient \insertCite{hill1968TAG}{GUSLD} defined as
#' \eqn{r^2=D^2/\sqrt{p_{A_1}p_{A_2}(1-p_{A_1})(1-p_{A_2})}}}
#' }
#' where \eqn{p_{A_1A_2}}  is the probability of observing a haplotype containing the reference allele at both loci,
#' \eqn{p_{A_1}} is the major allele frequency at the first SNP and \eqn{p_{A_2}} is the allele
#' frequency at the second SNP.
#' One can also specify their own LD measure (see examples below).
#'
#' @param pA1 Allele frequency for the first SNP
#' @param pA2 Allele frequency for the second SNP
#' @param D Disequilibrium coefficient
#'
#' @return A numeric value giving the value of the LD measure
#' @name LDmeasure
#' @author Timothy P. Bilton
#' @references
#' \insertAllCited{}
#' @examples
#' Dcoef(0.25,0.25,0.1)
#' Dprime(0.25,0.25,0.1)
#' r2(0.25,0.25,0.1)
#'
#' ## One can also define there own LD measure which GUSLD can use
#' ## e.g., the correlation coefficient
#' LD_r <- function(pA1,pA2,D){
#' return( D/sqrt((prod(c(pA1,pA2,1-c(pA1,pA2))))) )
#' }

#' @export Dcoef
#' @rdname LDmeasure
Dcoef <- function(pA1,pA2,D){
  return(D)
}

#' @export Dprime
#' @rdname LDmeasure
Dprime <- function(pA1,pA2,D){
  C1 = max(-prod(c(pA1,pA2)),-prod(1-c(pA1,pA2))); C2 = min((1-pA1)*pA2,pA1*(1-pA2))
  return(D/ifelse(D<0,C1,C2))
}

#' @export r2
#' @rdname LDmeasure
r2 <- function(pA1,pA2,D){
  return(D^2/(prod(c(pA1,pA2,1-c(pA1,pA2)))))
}


