#' @export GUSLD


GUSLD <- function(R6obj, set1=NULL, set2=NULL, nClust=3, LDmeasure, toFile=F, resfile=NULL, file="LDresults"){

  ## do some checks
  if(!all(class(R6obj) %in% c("UR","RA","R6")))
    stop("First argumented supplied is not of class 'UR', 'R6' and 'RA'")
  if(!is.vector(LDmeasure) || !is.character(LDmeasure) || length(LDmeasure) == 0)
    stop("No LD measure has been specified")
  if(!is.vector(toFile) || length(toFile) != 1 || !is.logical(toFile))
    stop("Argument for writing results is invalid")
  if(toFile){
    file = file.path(file,".txt")
    if(file.exists(file))
      stop("Output file already exists. Please specific a different file name")
  }
  nSnps <- R6obj$.__enclos_env__$private$nSnps
  if(is.null(set1))
    set1 <- 1:nSnps
  else if(!is.numeric(set1) || !is.vector(set1) || round(set1) == set1 ||
          any(set1 < 0) || any(set1 > nSnps))
    stop("Argument for indexing of first SNP set is invalid")
  else
    set1 <- sort(unique(set1))
  if(is.null(set2))
    set2 <- 1:nSnps
  else if(!is.numeric(set2) || !is.vector(set2) || round(set2) == set2 ||
          any(set2 < 0) || any(set2 > nSnps))
    stop("Argument for indexing of second SNP set is invalid")
  else
    set2 <- sort(unique(set2))
  nSnps1 <- length(set1)
  nSnps2 <- length(set2)
  ## estimate the pairwise LD
  res <- foreach(snp1=iter(1:nSnps1),.combine=comb,.multicombine=TRUE) %dopar% {
    LDvec <- c(replicate(length(LDmeasure),numeric(nSnps2),simplify=F))
    for(snp2 in seq_len(snp1-1)){
      ind <- c(snp1,snp2)
      MLE <- try(optim(par = 0, fn = "ll_gusld", LD=0,
                   p=R6obj$.__enclos_env__$private$p[ind],
                   ep=R6obj$.__enclos_env__$private$ep[ind],
                   ref=R6obj$.__enclos_env__$private$ref[,ind],
                   alt=R6obj$.__enclos_env__$private$alt[,ind],
                   nInd=R6obj$.__enclos_env__$private$nInd))
      ## check that the estimation process worked
      if(class(MLE)=="try-error")
        MLE = NA
      else
        MLE = MLE$par
      ## generate the estimates
      depth <- R6obj$.__enclos_env__$private$ref[,ind] + R6obj$.__enclos_env__$private$alt[,ind]
      N <- apply(depth, 1, function(x) all(!is.na(x)))
      D_hat <- (2*inv.logit(MLE[1])-1)*sqrt(pA1_hat*pA2_hat*(1-pA1_hat)*(1-pA2_hat))
      D_hat <- D_hat*(2*N/(2*N-1))
      C1hat <- max(-(pA1_hat*pA2_hat),-prod(1-c(pA1_hat,pA2_hat))); C2hat <- min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
      D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
      for(meas in 1:length(LDmeasure)) {
        LDvec[[meas]][snp2] <- get(LDmeasure[[meas]])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
      }
    }
    return(LDvec,round)
  }
  ## check what to do with the output
  if(toFile)
    write.table(res,file = file)
  else{
    return(res)
  }
}




}
