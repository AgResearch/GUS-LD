% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDmeasures.R
\name{LDmeasure}
\alias{LDmeasure}
\alias{Dcoef}
\alias{Dprime}
\alias{r2}
\title{Linkage Disequilibrium (LD) Measures}
\usage{
Dcoef(pA1, pA2, D)

Dprime(pA1, pA2, D)

r2(pA1, pA2, D)
}
\arguments{
\item{pA1}{Allele frequency for the first SNP}

\item{pA2}{Allele frequency for the second SNP}

\item{D}{Disequilibrium coefficient}
}
\value{
A numeric value giving the value of the LD measure
}
\description{
Functions for computing LD measures for predefined (or estimated) allele frequency values and
a predefined (or estimated) disequilibrium coefficient value.
}
\details{
The LD measures available in GUSLD are:
\describe{
\item{\code{Dcoef}}{The disequilibrium coefficient \insertCite{lewontin1960evolution}{GUSLD}
defined as \eqn{D=p_{A_1A_2}-p_{A_1}p_{A_2}}}
\item{\code{Dprime}}{The normalized disequilibrium coefficient \insertCite{lewontin1964genetics}{GUSLD} defined as
\eqn{D'=D/D_{max}} where \eqn{D_{max} = min(p_{A_1}p_{A_2},(1-p_{A_1})(1-p_{A_2}))} if \eqn{D>0} and
\eqn{D_{max} = min(p_{A_1}(1-p_{A_2}),(1-p_{A_1})p_{A_2})}}
\item{\code{r2}}{The squared correlation coefficient \insertCite{hill1968TAG}{GUSLD} defined as
\eqn{r^2=D^2/\sqrt{p_{A_1}p_{A_2}(1-p_{A_1})(1-p_{A_2})}}}
}
where \eqn{p_{A_1A_2}}  is the probability of observing a haplotype containing the reference allele at both loci,
\eqn{p_{A_1}} is the major allele frequency at the first SNP and \eqn{p_{A_2}} is the allele
frequency at the second SNP.
One can also specify their own LD measure (see examples below).
}
\examples{
Dcoef(0.25,0.25,0.1)
Dprime(0.25,0.25,0.1)
r2(0.25,0.25,0.1)

## One can also define there own LD measure which GUSLD can use
## e.g., the correlation coefficient
LD_r <- function(pA1,pA2,D){
return( D/sqrt((prod(c(pA1,pA2,1-c(pA1,pA2))))) )
}
}
\references{
\insertAllCited{}
}
\author{
Timothy P. Bilton
}
