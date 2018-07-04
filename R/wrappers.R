#' @useDynLib GUSLD



## likelihood function for LD estimation
ll_gusld <- function(LD, p, ep, ref, alt, nInd){
  epMat <- matrix(ep,nrow=nInd, ncol=2)
  AA <- (1-epMat)^ref*epMat^alt
  AB <- (1/2)^(ref+alt)
  BB <- epMat^ref*(1-epMat)^alt
  ll <- .Call("ll_gusld_c",LD, p, ep, AA, AB, BB, nInd)
  return(-ll)
}
