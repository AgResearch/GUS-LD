##########################################################################
# Genotyping Uncertainty with Sequencing data - Linkage Disequilibrium (GUSLD)
# Copyright 2017-2025 AgResearch Ltd. <timothy.bilton@agresearch.co.nz>
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

#' Compute pairwise LD estimates using GUSLD
#'
#' Function computes the LD estimates for high-throughput sequencing data using the methodology
#' by \insertCite{bilton2018genetics2;textual}{GUSbase}.
#'
#' Argument \code{SNPpairs} must be an integer matrix with two columns, where the first column specify the
#' indices of the SNPs and the rows correspond to each SNP pair for which LD is to be computed.
#' If \code{SNPpairs=NULL}, then LD is estimated for all the SNP paris in the
#' dataset.
#'
#' To reduce computation time, the calculations are performed in parallel using the \code{\link[foreach]{foreach}}
#' function. The \code{nClust} argument specifies the number of cores to use in the parallelization. Note: do not
#' set this argument to more than the number of cores available (or bad things can happen!).
#'
#' @param URobj An object of class UR created by the \code{\link{makeUR}} function.
#' @param SNPpairs Matrix specifying the pairs of SNPs in which to calculate pairwise LD. See below for details.
#' @param indsubset Vector of integer indices corresponding to individuals to be used in the LD calculations.
#' Used to specify a subset of individuals of the data for which LD is calculated.
#' @param nClust Integer number specifying the number of cores to use in the paralellization.
#' @param LDmeasure Character vector specifying which LD measures to calculate. Predefined
#' measures are \code{\link{Dcoef}}, \code{\link{Dprime}} and \code{\link{r2}} (see details below).
#' One can also specify their own LD measure (see examples).
#' @param filename Character string giving the name of the file to write the LD results to. If NULL,
#' the function returns the LD results instead of writing them to a file.
#' @param dp Integer number specifying the number of decimal places to round the LD results.
#' Note: only works when \code{file} is not NULL.
#'
#' @return
#' A matrix object with named columns is returned where the columns are:
#' \itemize{
#' \item [\emph{LDmeasure}]: The LD estmate for the measure [\emph{LDmeasure}]. Note: if multiple LD measures have been specified,
#' then there will be a column for each LD measure
#' \item CHROM_SNP1: The chromosome number of the first SNP
#' \item CHROM_SNP2: The chromosome number of the second SNP
#' \item POS_SNP1: The position (in base pairs) of the first SNP on the chromosome
#' \item POS_SNP2: The position (in base pairs) of the second SNP on the chromosome
#' \item FREQ_SNP1: The allele frequency estimate of the first SNP (as used in the calculation of the LD measure)
#' \item FREQ_SNP2: The allele frequency estimate of the second SNP (as used in the calculation of the LD measure)
#' \item ERR_SNP1: The sequencing error estimate of the first SNP (as used in the calculation of the LD measure)
#' \item ERR_SNP2: The sequencing error estimate of the second SNP (as used in the calculation of the LD measure)
#' }
#' If \code{filename=NULL} then the matrix object is returned to the R workspace, otherwise if \code{filename} is specified
#' and corresponds to a vaild name, then the matrix object is written to a file called \emph{filename_GUSMap.txt}.
#' @author Timothy P. Bilton
#' @references
#' \insertRef{bilton2018genetics2}{GUSbase}
#' @examples
#'
#' ## Read in the deer data that accompanies the package
#' deerfile <- deerData()
#' rafile <- VCFtoRA(deerfile)
#' deer <- readRA(rafile)
#'
#' ## Create unrelated population
#' ur <- makeUR(deer)
#'
#' ###### LD estimation #######
#' ## Estimate all the pairwise LD
#' LDres <- GUSLD(ur)
#' head(LDres)
#'
#' ## write results to file
#' GUSLD(ur, file="results")
#'
#' ## Specifying SNP pairs explicitly
#' # block vs block
#' pairs <- as.matrix(expand.grid(1:10,11:20))
#' LDres <- GUSLD(ur, SNPpairs=pairs)
#' head(LDres)
#'
#' # Five SNPs before and after
#' temp <- rep(1:37,rep(5,37))
#' pairs <- cbind(temp, temp + 1:5)
#' pairs <- pairs[-which(pairs[,2] > 38),]
#' LDres <- GUSLD(ur, SNPpairs=pairs)
#' head(LDres)
#'
#' ##### Non-standard LD measure: #######
#' ## Define LD measure as the correlation coefficient
#' LD_r <- function(pA1,pA2,D){
#' return( D/sqrt((prod(c(pA1,pA2,1-c(pA1,pA2))))) )
#' }
#' LDres <- GUSLD(ur, LDmeasure="LD_r")
#' head(LDres)
#'
#' @export
GUSLD <- function(URobj, SNPpairs=NULL, indsubset=NULL, nClust=2, LDmeasure="r2",
                  filename=NULL, dp=4){

  ## do some checks
  if(!all(class(URobj) %in% c("UR","RA","R6")))
    stop("First argumented supplied is not of class 'UR', 'R6' and 'RA'")
  if(!is.vector(LDmeasure) | !is.character(LDmeasure) | length(LDmeasure) == 0)
    stop("Argument for which LD measures to use is invalid.")
  if(!is.numeric(nClust) | length(nClust) != 1 | round(nClust) != nClust | nClust < 1)
    stop("Argument for the number of cores in the parallel processing is invalid. Should be a positive integer number")
  else{
    for(measure in LDmeasure){
      if(!exists(measure) | !is.function(get(measure)))
        stop(paste("LD measure",measure,"is not defined"))
    }
  }
  if(is.null(indsubset))
    indsubset <- 1:URobj$.__enclos_env__$private$nInd
  else if(!(is.vector(indsubset) && is.numeric(indsubset) && all(round(indsubset) == indsubset) &&
            min(indsubset) > 0 && max(indsubset) < (URobj$.__enclos_env__$private$nInd+1))){
    stop(paste0("Argument for individuals is not a integer vector between 1 and the number fo individuals",
                URobj$.__enclos_env__$private$nInd))
  } else indsubset <- unique(indsubset)
  nind <- length(indsubset)
  ## Check if we write to filename
  writeFile = FALSE
  if(!is.null(filename)){
    if(!is.vector(filename) | !is.character(filename) | length(filename) != 1)
      stop("Specified file name is invalid")
    else {
      outfilename <- paste0(tail(strsplit(filename,split=.Platform$file.sep)[[1]],1),"_GUSLD.txt")
      writeFile = TRUE
      outpath <- GUSbase::dts(normalizePath("./", winslash=.Platform$file.sep, mustWork=T))
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
    #cl <- parallel::makeCluster(nClust)
    doParallel::registerDoParallel(nClust)
    nSnps <- length(snps)
    ## estimate the pairwise LD
    res <- foreach::foreach(snp1=1:nSnps,.combine=GUSbase:::comb_mat, .export=LDmeasure,
                   .multicombine=TRUE) %dopar% {
      LDvec <- c(replicate(length(LDmeasure) + 2,numeric(nSnps),simplify=F))
      for(snp2 in seq_len(snp1-1)){
        ind <- snps[c(snp1,snp2)]
        temp <- URobj$.__enclos_env__$private$pfreq[ind]
        pA1_hat = temp[1]
        pA2_hat = temp[2]
        C1hat = max(-prod(c(pA1_hat,pA2_hat)),-prod(1-c(pA1_hat,pA2_hat)))
        C2hat = min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
        ref_temp <- URobj$.__enclos_env__$private$ref[indsubset,ind]
        alt_temp <- URobj$.__enclos_env__$private$alt[indsubset,ind]
        depth <- ref_temp + alt_temp
        N <- sum(apply(depth, 1, function(x) all(x > 0)))
        Neff <- sum(1-(1/2)^(apply(depth[,1:2], 1, function(x) min(x))))
        LDvec[[length(LDmeasure) + 1]][snp2] <- N
        LDvec[[length(LDmeasure) + 2]][snp2] <- Neff
        if(N > 2){
          MLE <- try(stats::optimize(f = ll_gusld, tol=1e-20,
                           lower=C1hat, upper=C2hat, p=c(pA1_hat,pA2_hat),
                           ep=URobj$.__enclos_env__$private$ep[ind],
                           ref=URobj$.__enclos_env__$private$ref[indsubset,ind],
                           alt=URobj$.__enclos_env__$private$alt[indsubset,ind],
                           nInd=nind))
          ## check that the estimation process worked
          if(class(MLE)=="try-error")
            D_hat = NA
          else
            D_hat = MLE$minimum
          ## generate the estimates
          D_hat <- D_hat*(2*N/(2*N-1))
          D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
          for(meas in 1:length(LDmeasure)) {
            LDvec[[meas]][snp2] <- get(LDmeasure[meas])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
          }
        }
        else{
          for(meas in 1:length(LDmeasure)) {
            LDvec[[meas]][snp2] <- NA
          }
        }
      }
      return(LDvec)
    }
    #parallel::stopCluster(cl)
    ## format result to write to file
    # if(writeFile){
      indx <- which(upper.tri(res[[1]]), arr.ind=T)
      out <- as.data.frame(matrix(nrow=nrow(indx), ncol=length(LDmeasure)+10))
      #out[,1:2] <- indx
      for(meas in 1:length(LDmeasure))
        out[meas] <- format(round(res[[meas]][indx],dp),scientific = FALSE,drop0trailing = TRUE)
      out[length(LDmeasure)+1:4] <-
        c(URobj$.__enclos_env__$private$chrom[snps][indx],
          URobj$.__enclos_env__$private$pos[snps][indx])
      out[length(LDmeasure)+5:8] <- format(c(
          round(URobj$.__enclos_env__$private$pfreq[snps][indx],dp),
          round(URobj$.__enclos_env__$private$ep[snps][indx],dp)))
      out[length(LDmeasure) + 9] <- res[[length(LDmeasure) + 1]][indx]
      out[length(LDmeasure) + 10] <- res[[length(LDmeasure) + 2]][indx]
      colnames(out) <- c(LDmeasure, "CHROM_SNP1","CHROM_SNP2","POS_SNP1","POS_SNP2",
                         "FREQ_SNP1","FREQ_SNP2","ERR_SNP1","ERR_SNP2","N","Neff")
      if(writeFile){
        data.table::fwrite(out, file = outfilename, quote=FALSE, nThread = nClust)
        cat("Name of LD results file:    \t",outfilename,"\n")
        cat("Location of LD results file:\t",outpath,"/\n\n", sep = "")
      } else{
    # } else{
    #   for(meas in 1:length(res)){
    #     diag(res[[meas]]) <- get(LDmeasure[meas])(pA1=0.5,pA2=0.5,D=0.25)
    #     res[[meas]][lower.tri(res[[meas]])] <- t(res[[meas]])[lower.tri(res[[meas]])]
    #     #rownames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps]
    #     #colnames(res[[meas]]) <- URobj$.__enclos_env__$private$SNP_Names[snps]
    #   }
    #   res[[length(res)+1]] <- URobj$.__enclos_env__$private$pfreq[snps]
    #   res[[length(res)+1]] <- URobj$.__enclos_env__$private$ep[snps]
    #   res[[length(res)+1]] <- URobj$.__enclos_env__$private$SNP_Names[snps]
    #   names(res) <- c(LDmeasure,"p","ep","SNP_Names")
      return(out)
    }
  } else{ ## Case 2: SNP pairs specified
    ## check input for SNPpairs matrix
    if(!is.numeric(SNPpairs) | !is.matrix(SNPpairs) | ncol(SNPpairs) != 2 | round(SNPpairs) != SNPpairs |
       min(SNPpairs) < 1 | max(SNPpairs) > URobj$.__enclos_env__$private$nSnps)
      stop("Input for SNPpairs is invalid. Should a integer matix with 2 columns.")
    # check for duplicates
    SNPpairs <- unique(SNPpairs, MARGIN=1)
    ## remove pairs where the SNP is the same
    removeSNPs <- which(SNPpairs[,1] == SNPpairs[,2])
    if(length(removeSNPs) > 1)
      SNPpairs <- SNPpairs[-removeSNPs,]
    ## Work out the number of pairs
    npairs <- nrow(SNPpairs)
    ## Set up the clusters
    #cl <- parallel::makeCluster(nClust)
    doParallel::registerDoParallel(nClust)
    ## compute the LD for specified SNP pairs
    res <- foreach::foreach(pair=1:npairs, .combine="rbind", .export=LDmeasure) %dopar% {
      LDvec <- rep(NA,length(LDmeasure))
      ind <- SNPpairs[pair,]
      temp <- URobj$.__enclos_env__$private$pfreq[ind]
      pA1_hat = temp[1]
      pA2_hat = temp[2]
      C1hat = max(-prod(c(pA1_hat,pA2_hat)),-prod(1-c(pA1_hat,pA2_hat)))
      C2hat = min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
      ref_temp <- URobj$.__enclos_env__$private$ref[indsubset,ind]
      alt_temp <- URobj$.__enclos_env__$private$alt[indsubset,ind]
      depth <- ref_temp + alt_temp
      N <- sum(apply(depth, 1, function(x) all(x > 0)))
      Neff <- sum(1-(1/2)^(apply(depth[,1:2], 1, function(x) min(x))))
      if(N > 2){
        MLE <- try(stats::optimize(f = ll_gusld, tol=1e-20,
                           lower=C1hat, upper=C2hat, p=c(pA1_hat,pA2_hat),
                           ep=URobj$.__enclos_env__$private$ep[ind],
                           ref=ref_temp, alt=alt_temp,
                           nInd=nind))
        ## check that the estimation process worked
        if(class(MLE)!="try-error"){
          D_hat = MLE$minimum
          ## generate the estimates
          D_hat <- D_hat*(2*N/(2*N-1))
          D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))
          for(meas in 1:length(LDmeasure)){
            LDvec[meas] <- get(LDmeasure[meas])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
          }
        }
      } else{
        for(meas in 1:length(LDmeasure)){
          LDvec[meas] <- NA
        }
      }
      return(c(LDvec,N,Neff))
    }
    #parallel::stopCluster(cl)
    res <- unname(res)
    out <- as.data.frame(matrix(nrow=npairs, ncol=length(LDmeasure) + 10))
    out[1:length(LDmeasure)] <- format(round(res[,1:length(LDmeasure)],dp), scientific = FALSE,drop0trailing = TRUE)
    out[length(LDmeasure) + 1:4] <- c(URobj$.__enclos_env__$private$chrom[SNPpairs],
                                      URobj$.__enclos_env__$private$pos[SNPpairs])
    out[length(LDmeasure) + 5:8] <- format(c(
      round(URobj$.__enclos_env__$private$pfreq[SNPpairs],dp),
      round(URobj$.__enclos_env__$private$ep[SNPpairs],dp)
    ),scientific = FALSE,drop0trailing = TRUE)
    out[length(LDmeasure) + 9] <- res[,length(LDmeasure) + 1]
    out[length(LDmeasure) + 10] <- res[,length(LDmeasure) + 2]
    names(out) <- c(LDmeasure, "CHROM_SNP1","CHROM_SNP2","POS_SNP1","POS_SNP2",
                    "FREQ_SNP1","FREQ_SNP2","ERR_SNP1","ERR_SNP2","N","Neff")
    ## return the results
    if(writeFile)
      data.table::fwrite(out,file = outfilename, quote = FALSE, nThread = nClust)
    else
      return(out)
  }
}
