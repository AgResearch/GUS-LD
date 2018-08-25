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





