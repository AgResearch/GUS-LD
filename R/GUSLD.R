#' @export GUSLD
#' @examples
#' ## Define LD measure as the correlation coefficient
#' LD_r <- function(pA1,pA2,D){
#' return( D/sqrt((prod(c(pA1,pA2,1-c(pA1,pA2))))) )
#' }
#'
#'

GUSLD <- function(URobj, SNPpairs=NULL, indset=NULL, nClust=3, LDmeasure="r2",
                  file=NULL){

  ## do some checks
  if(!all(class(URobj) %in% c("UR","RA","R6")))
    stop("First argumented supplied is not of class 'UR', 'R6' and 'RA'")
  if(!is.vector(LDmeasure) || !is.character(LDmeasure) || length(LDmeasure) == 0)
    stop("Argument for which LD measures to use is invalid.")
  if(!is.numeric(nClust) || length(nClust) != 1 || round(nClust) != nClust || nClust < 1)
    stop("Argument for the number of cores in the parallel processing is invalid. Should be a positive integer number")
  else{
    for(measure in LDmeasure){
      if(!exists(measure) || !is.function(get(measure)))
        stop(paste("LD measure",measure,"is not defined"))
    }
  }
  if(is.null(indset))
    indset <- 1:URobj$.__enclos_env__$private$nInd
  else if(!(is.vector(indset) && is.numeric(indset) && all(round(indset) == indset) &&
            min(indset) > 0 && max(indset) < (URobj$.__enclos_env__$private$nInd+1))){
    stop(paste0("Argument for individuals is not a integer vector between 1 and the number fo individuals",
                URobj$.__enclos_env__$private$nInd))
  } else indset <- unique(indset)
  nind <- length(indset)
  ## Check if we write to file
  writeFile = FALSE
  if(!is.null(file)){
    if(!is.vector(file) || !is.character(file) || length(file) != 1)
      stop("Specified file name is invalid")
    else {
      file <- paste0("./",file,"_GUSLD.txt")
      writeFile = TRUE
      #if(file.access(file) == 0)
      #  writeFile = TRUE
      #else
      #  warning("Cannot write to specified location. Returning output to console.")
    }
  }
  ## Case 1: All pairs
  if(is.null(SNPpairs)){
    snps <- 1:URobj$.__enclos_env__$private$nSnps
    ## Set up the clusters
    cl <- makeCluster(nClust)
    registerDoSNOW(cl)
    nSnps <- length(snps)
    ## estimate the pairwise LD
    res <- foreach(snp1=iter(1:nSnps),.combine=GUSbase:::comb_mat,.multicombine=TRUE) %dopar% {
      LDvec <- c(replicate(length(LDmeasure),numeric(nSnps),simplify=F))
      for(snp2 in seq_len(snp1-1)){
        ind <- snps[c(snp1,snp2)]
        temp <- URobj$.__enclos_env__$private$pfreq[ind]
        pA1_hat = temp[1]
        pA2_hat = temp[2]
        C1hat = max(-prod(c(pA1_hat,pA2_hat)),-prod(1-c(pA1_hat,pA2_hat)))
        C2hat = min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
        MLE <- try(optim(par = 0, fn = ll_gusld, method="Brent",
                         lower=C1hat, upper=C2hat, p=c(pA1_hat,pA2_hat),
                         ep=URobj$.__enclos_env__$private$ep[ind],
                         ref=URobj$.__enclos_env__$private$ref[indset,ind],
                         alt=URobj$.__enclos_env__$private$alt[indset,ind],
                         nInd=nind))
        ## check that the estimation process worked
        if(class(MLE)=="try-error")
          MLE = NA
        else
          D_hat = MLE$par
        ## generate the estimates
        depth <- URobj$.__enclos_env__$private$ref[indset,ind] + URobj$.__enclos_env__$private$alt[indset,ind]
        N <- sum(apply(depth, 1, function(x) all(!is.na(x))))
        D_hat <- D_hat*(2*N/(2*N-1))
        D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
        for(meas in 1:length(LDmeasure)) {
          LDvec[[meas]][snp2] <- get(LDmeasure[meas])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
        }
      }
      return(LDvec)
    }
    stopCluster(cl)
    ## format result to write to file
    if(writeFile){
      indx <- which(upper.tri(res[[1]]), arr.ind=T)
      out <- as.data.frame(matrix(nrow=nrow(indx), ncol=length(LDmeasure)+8))
      #out[,1:2] <- indx
      for(meas in 1:length(LDmeasure))
        out[meas] <- round(res[[meas]][indx],4)
      out[length(LDmeasure)+1:4] <-
        c(URobj$.__enclos_env__$private$chrom[snps][indx],
          URobj$.__enclos_env__$private$pos[snps][indx])
      out[length(LDmeasure)+5:8] <- format(c(
          round(URobj$.__enclos_env__$private$pfreq[snps][indx],4),
          format(round(URobj$.__enclos_env__$private$ep[snps][indx],4))),scientific = FALSE,drop0trailing = TRUE)
      colnames(out) <- c(LDmeasure, "CHROM_SNP1","CHROM_SNP2","POS_SNP1","POS_SNP2",
                         "FREQ_SNP1","FREQ_SNP2","ERR_SNP1","ERR_SNP2")
      data.table::fwrite(out, file = file, quote=FALSE, nThread = nClust)
    } else{
      for(meas in 1:length(res)){
        diag(res[[meas]]) <- get(LDmeasure[meas])(pA1=0.5,pA2=0.5,D=0.25)
        res[[meas]][lower.tri(res[[meas]])] <- t(res[[meas]])[lower.tri(res[[meas]])]
        #rownames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps]
        #colnames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps]
      }
      res[[length(res)+1]] <- URobj$.__enclos_env__$private$pfreq[snps]
      res[[length(res)+1]] <- URobj$.__enclos_env__$private$ep[snps]
      res[[length(res)+1]] <- URobj$.__enclos_env__$private$SNP_Names[snps]
      res[[length(res)+1]] <- URobj$.__enclos_env__$private$SNP_Names[snps]
      names(res) <- c(LDmeasure,"p","ep","SNP Names (rows)","SNP Names (columns)")
      return(res)
    }
  } else{ ## Case 2: SNP pairs specified
    ## check input for SNPpairs matrix
    if(!is.numeric(SNPpairs) || !is.matrix(SNPpairs) || ncol(SNPpairs) != 2 || round(SNPpairs) != SNPpairs ||
       min(SNPpairs) < 1 || max(SNPpairs) > URobj$.__enclos_env__$private$nSnps)
      stop("Input for SNPpairs is invalid. Should a integer matix with 2 columns.")
    # check for duplicates
    SNPpairs <- unique(SNPpairs, MARGIN=1)
    ## Work out the number of pairs
    npairs <- nrow(SNPpairs)
    ## Set up the clusters
    cl <- makeCluster(nClust)
    registerDoSNOW(cl)
    ## compute the LD for specified SNP pairs
    res <- foreach(pair=iter(1:npairs), .combine="rbind") %dopar% {
      LDvec <- rep(NA,length(LDmeasure))
      ind <- SNPpairs[pair,]
      temp <- URobj$.__enclos_env__$private$pfreq[ind]
      pA1_hat = temp[1]
      pA2_hat = temp[2]
      C1hat = max(-prod(c(pA1_hat,pA2_hat)),-prod(1-c(pA1_hat,pA2_hat)))
      C2hat = min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
      MLE <- try(optim(par = 0, fn = ll_gusld, method="Brent",
                         lower=C1hat, upper=C2hat, p=c(pA1_hat,pA2_hat),
                         ep=URobj$.__enclos_env__$private$ep[ind],
                         ref=URobj$.__enclos_env__$private$ref[indset,ind],
                         alt=URobj$.__enclos_env__$private$alt[indset,ind],
                         nInd=nind))
      ## check that the estimation process worked
      if(class(MLE)!="try-error"){
        D_hat = MLE$par
        ## generate the estimates
        depth <- URobj$.__enclos_env__$private$ref[indset,ind] + URobj$.__enclos_env__$private$alt[indset,ind]
        N <- sum(apply(depth, 1, function(x) all(!is.na(x))))
        D_hat <- D_hat*(2*N/(2*N-1))
        D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
        for(meas in 1:length(LDmeasure)){
          LDvec[meas] <- get(LDmeasure[meas])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
        }
      }
      return(LDvec)
    }
    stopCluster(cl)
    out <- as.data.frame(matrix(nrow=npairs, ncol=length(LDmeasure) + 8))
    out[1:length(LDmeasure)] <- round(res,4)
    out[length(LDmeasure) + 1:4] <- c(URobj$.__enclos_env__$private$chrom[SNPpairs],
                                      URobj$.__enclos_env__$private$pos[SNPpairs])
    out[length(LDmeasure) + 5:8] <- format(c(
      round(URobj$.__enclos_env__$private$pfreq[SNPpairs],4),
      round(URobj$.__enclos_env__$private$ep[SNPpairs],4)
    ),scientific = FALSE,drop0trailing = TRUE)
    names(out) <- c(LDmeasure, "CHROM_SNP1","CHROM_SNP2","POS_SNP1","POS_SNP2",
                    "FREQ_SNP1","FREQ_SNP2","ERR_SNP1","ERR_SNP2")
    ## return the results
    if(writeFile)
      data.table::fwrite(out,file = file, quote = FALSE, nThread = nClust)
    else
      return(out)
  }
}


#   if(!is.null(file)){
#       if(is.null(SNPpairs) && (is.null(SNPsets) || (is.list(SNPsets) && length(SNPsets) == 1))
#       ## Symmetric matrix for all pairs
#       if(is.null(SNPsets) || (is.list(SNPsets) && length(SNPsets) == 1 &&
#                               is.vector(SNPsets[[1]]) && round(SNPsets[[1]]) == SNPsets[[1]])){
#         if(!is.null(SNPsets)){
#           snps <- sort(SNPsets[[1]])
#           if(any(snps) < 1 || any(snps) > URobj$.__enclos_env__$private$nSnps || any(duplicated(snps)))
#             stop("SNPset vector is invalid")
#         }
#         else
#           snps <- 1:URobj$.__enclos_env__$private$nSnps
#
#         ## Set up the clusters
#         cl <- makeCluster(nClust)
#         registerDoSNOW(cl)
#         nSnps <- length(snps)
#         ## estimate the pairwise LD
#         res <- foreach(snp1=iter(1:nSnps),.combine=GUSbase:::comb_mat,.multicombine=TRUE) %dopar% {
#           LDvec <- c(replicate(length(LDmeasure),numeric(nSnps),simplify=F))
#           for(snp2 in seq_len(snp1-1)){
#             ind <- snps[c(snp1,snp2)]
#             temp <- URobj$.__enclos_env__$private$pfreq[ind]
#             pA1_hat = temp[1]
#             pA2_hat = temp[2]
#             C1hat = max(-prod(c(pA1_hat,pA2_hat)),-prod(1-c(pA1_hat,pA2_hat)))
#             C2hat = min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
#             MLE <- try(optim(par = 0, fn = ll_gusld, method="Brent",
#                              lower=C1hat, upper=C2hat, p=c(pA1_hat,pA2_hat),
#                              ep=URobj$.__enclos_env__$private$ep[ind],
#                              ref=URobj$.__enclos_env__$private$ref[,ind],
#                              alt=URobj$.__enclos_env__$private$alt[,ind],
#                              nInd=URobj$.__enclos_env__$private$nInd))
#             ## check that the estimation process worked
#             if(class(MLE)=="try-error")
#               MLE = NA
#             else
#               D_hat = MLE$par
#             ## generate the estimates
#             depth <- URobj$.__enclos_env__$private$ref[,ind] + URobj$.__enclos_env__$private$alt[,ind]
#             N <- sum(apply(depth, 1, function(x) all(!is.na(x))))
#             D_hat <- D_hat*(2*N/(2*N-1))
#             D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
#             for(meas in 1:length(LDmeasure)) {
#               LDvec[[meas]][snp2] <- get(LDmeasure[[meas]])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
#             }
#           }
#           return(LDvec)
#         }
#         stopCluster(cl)
#         for(meas in 1:length(res)){
#           diag(res[[meas]]) <- 0
#           res[[meas]][upper.tri(res[[meas]])] <- t(res[[meas]])[upper.tri(res[[meas]])]
#           #rownames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps]
#           #colnames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps]
#         }
#         res[[length(res)+1]] <- URobj$.__enclos_env__$private$pfreq[snps]
#         res[[length(res)+1]] <- URobj$.__enclos_env__$private$ep[snps]
#         res[[length(res)+1]] <- URobj$.__enclos_env__$private$SNP_Names[snps]
#         res[[length(res)+1]] <- URobj$.__enclos_env__$private$SNP_Names[snps]
#         names(res) <- c(LDmeasure,"p","ep","SNP Names (rows)","SNP Names (columns)")
#         if(!is.null(file))
#           stop("to be implemented")
#         else
#           return(res)
#       }
#       ## Block matrix for estimation between two independent pair of SNPs
#       else if(is.list(SNPsets) && length(SNPsets) == 2 &&
#               is.vector(SNPsets[[1]]) && round(SNPsets[[1]]) == SNPsets[[1]] &&
#               is.vector(SNPsets[[2]]) && round(SNPsets[[2]]) == SNPsets[[2]]){
#         snps1 <- sort(SNPsets[[1]])
#         snps2 <- sort(SNPsets[[2]])
#         if(any(snps1) < 1 || any(snps1) > URobj$.__enclos_env__$private$nSnps || any(duplicated(snps1)))
#           stop("SNPset vector is invalid")
#         if(any(snps2) < 1 || any(snps2) > URobj$.__enclos_env__$private$nSnps || any(duplicated(snps2)))
#           stop("SNPset vector is invalid")
#         ## Set up the clusters
#         cl <- makeCluster(nClust)
#         registerDoSNOW(cl)
#         ## estimate the pairwise LD
#         res <- foreach(snp1=iter(snps1),.combine=GUSbase:::comb_mat,.multicombine=TRUE) %dopar% {
#           LDvec <- c(replicate(length(LDmeasure),numeric(length(snps2)),simplify=F))
#           for(snp2 in snps2){
#             ind <- c(snp1,snp2)
#             temp <- URobj$.__enclos_env__$private$pfreq[ind]
#             pA1_hat = temp[1]
#             pA2_hat = temp[2]
#             C1hat = max(-prod(c(pA1_hat,pA2_hat)),-prod(1-c(pA1_hat,pA2_hat)))
#             C2hat = min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
#             MLE <- try(optim(par = 0, fn = ll_gusld, method="Brent",
#                              lower=C1hat, upper=C2hat, p=c(pA1_hat,pA2_hat),
#                              ep=URobj$.__enclos_env__$private$ep[ind],
#                              ref=URobj$.__enclos_env__$private$ref[,ind],
#                              alt=URobj$.__enclos_env__$private$alt[,ind],
#                              nInd=URobj$.__enclos_env__$private$nInd))
#             ## check that the estimation process worked
#             if(class(MLE)=="try-error")
#               MLE = NA
#             else
#               D_hat = MLE$par
#             ## generate the estimates
#             depth <- URobj$.__enclos_env__$private$ref[,ind] + URobj$.__enclos_env__$private$alt[,ind]
#             N <- sum(apply(depth, 1, function(x) all(!is.na(x))))
#             D_hat <- D_hat*(2*N/(2*N-1))
#             D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
#             for(meas in 1:length(LDmeasure)) {
#               LDvec[[meas]][which(snp2==snps2)] <- get(LDmeasure[[meas]])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
#             }
#           }
#           return(LDvec)
#         }
#         stopCluster(cl)
#         ## add the SNP names to the matrix
#         for(meas in 1:length(res)){
#           rownames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps1]
#           colnames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps2]
#         }
#         res[[length(res)+1]] <- list(Rows=URobj$.__enclos_env__$private$pfreq[snps1],
#                                       Columns=URobj$.__enclos_env__$private$pfreq[snps2])
#         res[[length(res)+1]] <- list(Rows=URobj$.__enclos_env__$private$ep[snps1],
#                                       Columns=URobj$.__enclos_env__$private$ep[snps2])
#         names(res) <- c(LDmeasure,"p","ep")
#         return(res)
#       }
#       else
#         stop("Invalid input for 'SNPsets' argument. Should be a list of length 2.")
#     }
#   }
#   ## write the results to a file
#   else if(output=="file"){
#     ## check file doesn't exist
#     if(output == "file"){
#       file = file.path(file,".txt")
#       if(file.exists(file))
#         stop("Output file already exists. Please specify a different file name")
#     }
#     if(is.null(SNPpairs)){
#       ## symmetic matrix between all paris of SNPs
#       if(is.null(SNPsets) || is.list(SNPsets) && length(SNPsets) == 1 &&
#          is.vector(SNPsets[[1]]) && round(SNPsets[[1]]) == SNPsets[[1]]){
#         if(is.null(SNPsets))
#           snps <- 1:URobj$.__enclos_env__$private$nSnps
#         else{
#           snps <- sort(SNPsets[[1]])
#           if(any(snps) < 1 || any(snps) > URobj$.__enclos_env__$private$nSnps || any(duplicated(snps)))
#             stop("SNPset vector is invalid")
#         }
#         # Every pairwise calculation
#         ## Set up the clusters
#         cl <- makeCluster(nClust)
#         registerDoSNOW(cl)
#         nSnps <- length(snps)
#         ## estimate the pairwise LD
#         res <- foreach(snp1=iter(1:nSnps),.combine=GUSbase:::comb_vec,.multicombine=TRUE) %dopar% {
#           LDvec <- c(replicate(length(LDmeasure)+1,numeric(length(seq_len(snp1-1))),simplify=F))
#           for(snp2 in seq_len(snp1-1)){
#             ind <- c(snp1,snp2)
#             temp <- URobj$.__enclos_env__$private$pfreq[ind]
#             pA1_hat = temp[1]
#             pA2_hat = temp[2]
#             C1hat = max(-prod(c(pA1_hat,pA2_hat)),-prod(1-c(pA1_hat,pA2_hat)))
#             C2hat = min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
#             MLE <- try(optim(par = 0, fn = ll_gusld, method="Brent",
#                              lower=C1hat, upper=C2hat, p=c(pA1_hat,pA2_hat),
#                              ep=URobj$.__enclos_env__$private$ep[ind],
#                              ref=URobj$.__enclos_env__$private$ref[,ind],
#                              alt=URobj$.__enclos_env__$private$alt[,ind],
#                              nInd=URobj$.__enclos_env__$private$nInd))
#             ## check that the estimation process worked
#             if(class(MLE)=="try-error")
#               MLE = NA
#             else
#               D_hat = MLE$par
#             ## generate the estimates
#             depth <- URobj$.__enclos_env__$private$ref[,ind] + URobj$.__enclos_env__$private$alt[,ind]
#             N <- sum(apply(depth, 1, function(x) all(!is.na(x))))
#             D_hat <- D_hat*(2*N/(2*N-1))
#             D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
#             for(meas in 1:length(LDmeasure)) {
#               LDvec[[meas]][snp2] <- get(LDmeasure[[meas]])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
#             }
#             LDvec[[meas+1]][snp2] <- c(URobj$.__enclos_env__$private$SNP_Names[ind[1]], pA1_hat, URobj$.__enclos_env__$private$ep[ind[1]],
#                                        URobj$.__enclos_env__$private$SNP_Names[ind[1]], pA2_hat, URobj$.__enclos_env__$private$ep[ind[2]])
#
#           }
#           return(LDvec)
#         }
#         stopCluster(cl)
#         return(do.call("cbind",res))
#       }
#       # Estimation between independent SNP sets
#       if(is.list(SNPsets) && length(SNPsets) == 2 && is.vector(SNPsets[[1]]) &&
#         round(SNPsets[[1]]) == SNPsets[[1]] && is.vector(SNPsets[[2]]) &&
#         round(SNPsets[[2]]) == SNPsets[[2]] && all(SNPsets[[1]] != SNPsets[[2]])){
#       }
#       else
#         stop("Invalid input for the 'SNPsets' argument. Should be a list of length 2")
#     }
#     ## Estimation for specified set of SNP pairs
#     else if(!is.null(SNPpairs) && is.null(SNPsets)){
#       if(is.matrix(SNPpairs) && is.numeric(SNPpairs) && all(round(SNPpairs)==SNPpairs) && nrow(SNPpairs) == 2){
#         pos_temp1 <- URobj$.__enclos_env__$private$SNP_Names[SNPpairs[,1]]
#         pos_temp2 <- URobj$.__enclos_env__$private$SNP_Names[SNPpairs[,2]]
#         SNPpairs <- cbind(SNPpairs,pos_temp1,pos_temp2)
#         rm(pos_temp1)
#         rm(pos_temp2)
#       }
#       else
#         stop("Invalid input for argument 'SNPpairs'. Should be an integer matrix with two columns")
#     }
#     else
#       stop("Something weird has happened. Should not get to this point.")
#   }
#   else
#     stop("Something weird has happened. Should not get to this point.")
# }
