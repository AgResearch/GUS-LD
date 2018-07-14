#' @export GUSLD


GUSLD <- function(R6obj, SNPpairs=NULL, SNPsets=NULL, nClust=3, LDmeasure=c("r2"), toFile=F, resfile=NULL, file="LDresults"){

  ## do some checks
  if(!all(class(R6obj) %in% c("UR","RA","R6")))
    stop("First argumented supplied is not of class 'UR', 'R6' and 'RA'")
  if(!is.vector(LDmeasure) || !is.character(LDmeasure) || length(LDmeasure) == 0)
    stop("No LD measure has been specified")
  else{
    for(measure in LDmeasure){
      if(!exists(measure) || !is.function(get(measure)))
        stop(paste("LD measure",measure,"is not defined")
    }
  }
  if(!is.vector(toFile) || length(toFile) != 1 || !is.logical(toFile))
    stop("Argument for writing results is invalid")
  if(toFile){
    file = file.path(file,".txt")
    if(file.exists(file))
      stop("Output file already exists. Please specific a different file name")
  }
  if(is.null(SNPpairs) && is.null(SNPsets)){ ## do every pairwise calculation
    ## Set up the clusters
    cl <- makeCluster(nClust)
    registerDoSNOW(cl)
    nSnps <- R6obj$.__enclos_env__$private$nSnps
    ## estimate the pairwise LD
    res <- foreach(snp1=iter(1:nSnps),.combine=comb,.multicombine=TRUE) %dopar% {
      LDvec <- c(replicate(length(LDmeasure),numeric(nSnps),simplify=F))
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
      return(LDvec)
    }
    stopCluster(cl)
    ## check what to do with the output
    if(toFile)
      write.table(res,file = file)
    else{
      return(res)
    }
  }
  else if(is.null(SNPpairs) && is.list(SNPsets) && length(SNPsets) == 2 &&
          is.vector(SNPsets[[1]]) && round(SNPsets[[1]]) == SNPsets[[1]] &&
          is.vector(SNPsets[[2]]) && round(SNPsets[[2]]) == SNPsets[[2]]){
    ## Set up the clusters
    cl <- makeCluster(nClust)
    registerDoSNOW(cl)
    nSnps1 <- length(SNPsets[[1]])
    nSnps2 <- length(SNPsets[[2]])
    ## estimate the pairwise LD
    res <- foreach(snp1=iter(1:nSnps1),.combine=comb,.multicombine=TRUE) %dopar% {
      LDvec <- c(replicate(length(LDmeasure),numeric(nSnps2),simplify=F))
      for(snp2 in nSnps2){
        ind <- c(SNPsets[snp1],SNPsets[snp2])
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
      return(LDvec)
    }
    stopCluster(cl)
    ## check what to do with the output
    if(toFile)
      write.table(res,file = file)
    else{
      return(res)
    }
  }
  else if(is.null(SNPsets) && is.matrix(SNPpairs) && all(round(SNPpairs)==SNPpairs) && nrow(SNPpairs) == 2){
    nSnps <- ncol(SNPpairs)
    ## Set up the clusters
    cl <- makeCluster(nClust)
    registerDoSNOW(cl)
    ## Do the calculations
    res <- foreach(snp=iter(1:nSnps),.combine=comb,.multicombine=TRUE) %dopar% {
      ind <- SNPpairs[snp,]
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
      LDvec <- c(ind,numeric(length(LDmeasure)))
      for(meas in 1:length(LDmeasure)) {
        LDvec[meas+2] <- get(LDmeasure[[meas]])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
      }
      return(LDvec)
    }
    stopCluster(cl)
    colnames(res) <- c("SNP1", "SNP2", LDmeasure)
    ## check what to do with the output
    if(toFile){
        write.table(res,file = file)
      else{
        return(res)
      }
    }
    else
      stop("Input for SNPpairs or/and SNPsets may be invalid")
  }
}

