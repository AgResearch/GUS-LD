
vcf <- deerData()
rafile <- VCFtoRA(vcf)
radata <- readRA(rafile,gform = "reference")

deer <- makeUR(radata, ploid = 2, filter = list(MAF=0.01, MISS=0.9))

## test out GUSLD

LDmat <- GUSLD(deer,nClust = 3, file="test",sf=4)
LDmat <- GUSLD(deer,SNPpairs = matrix(c(1:5,2:6), ncol=2),
               nClust = 3, file="test_sub", sf=6, LDmeasure = c("r2","LD_r"))


## Simulation

estMLEadjustedLik2 <- function(depth_Ref,depth_Alt,initialV){

  pAB <- choose(depth_Ref+depth_Alt,depth_Ref) * (1/2)^(depth_Ref + depth_Alt)
  N <- length(apply(depth_Ref==0&depth_Alt==0,1,any))

  ## check to see if the initial value starting values are not too close to zero
  if(initialV[1] > 0.999)
    initialV[1] = 0.999
  if(initialV[2] > 0.999)
    initialV[2] = 0.999

  ## likelihood function

  ll_wrapper <- function(logit_x, depth_Ref, depth_Alt, pAB){
    pA1=inv.logit(logit_x[2])
    pA2=inv.logit(logit_x[3])
    ep =inv.logit(logit_x[4])
    r = 2*inv.logit(logit_x[1]) - 1
    D = r*sqrt(pA1*pA2*(1-pA1)*(1-pA2))
    return(ll(c(pA1,pA2,D,ep), depth_Ref, depth_Alt, pAB))
  }
  ll <- function(x, depth_Ref, depth_Alt, pAB){
    pA1=x[1]; pA2=x[2]; D=x[3]
    pB1 = 1-pA1; pB2 = 1-pA2
    C1 <- max(-pA1*pA2,-pB1*pB2); C2 <- min(pB1*pA2,pA1*pB2)
    if(!((C1 <= D) & (D <= C2))) return(NA)
    else{
      ep = x[4]
      ## compute the probability matrices
      pAA <- choose(depth_Ref+depth_Alt,depth_Ref) * (1-ep)^depth_Ref * ep^depth_Alt
      pBB <- choose(depth_Ref+depth_Alt,depth_Ref) * (ep)^depth_Ref * (1-ep)^depth_Alt

      h1 <- pA1*pA2+D; h2 <- pA1*pB2-D; h3 <- pB1*pA2-D; h4 <- pB1*pB2+D
      p11 <- h1^2
      p21 <- 2*h1*h3
      p31 <- h3^2
      p12 <- 2*h1*h2
      p32 <- 2*h3*h4
      p13 <- h2^2
      p23 <- 2*h2*h4
      p33 <- h4^2
      p22 <- 2*h1*h4+2*h2*h3

      return(-(sum( log(
        pAA[,1]*pAA[,2]*p11 + pAA[,1]*pAB[,2]*p12 + pAA[,1]*pBB[,2]*p13 +
          pAB[,1]*pAA[,2]*p21 + pAB[,1]*pAB[,2]*p22 + pAB[,1]*pBB[,2]*p23 +
          pBB[,1]*pAA[,2]*p31 + pBB[,1]*pAB[,2]*p32 + pBB[,1]*pBB[,2]*p33
      )))  )
    }
  }
  ## Compute MLE
  MLE <- optim(logit(c(0.5,initialV[1:2],0.01)),ll_wrapper,
               method="Nelder-Mead", control=list(reltol=1e-20, maxit=2000), depth_Ref=depth_Ref, depth_Alt=depth_Alt, pAB=pAB)
  MLE2 <- optim(logit(c(0.5,initialV[1:2],0.0001)),ll_wrapper,
                method="Nelder-Mead", control=list(reltol=1e-20, maxit=2000), depth_Ref=depth_Ref, depth_Alt=depth_Alt, pAB=pAB)
  if(MLE$value < MLE2$value)
    MLE <- MLE$par
  else
    MLE <- MLE2$par
  pA1_hat = inv.logit(MLE[2])
  pA2_hat = inv.logit(MLE[3])
  D_hat <- (2*inv.logit(MLE[1])-1)*sqrt(pA1_hat*pA2_hat*(1-pA1_hat)*(1-pA2_hat))
  D_hat <- D_hat*(2*N/(2*N-1))
  ep_hat = inv.logit(MLE[4])
  ## Check that D is within its range
  C1hat <- max(-(pA1_hat*pA2_hat),-prod(1-c(pA1_hat,pA2_hat))); C2hat <- min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
  D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))

  ## Compute the estimates of Dprime and r2
  Dprime <- D_hat/ifelse(D_hat<0, abs(C1hat), C2hat)
  r2 <- D_hat^2/(pA1_hat*pA2_hat*(1-pA1_hat)*(1-pA2_hat))
  # return the estimates
  return(c(D_hat,Dprime,r2,pA1_hat,pA2_hat,ep_hat))
}

library(GUSLD)

nInd = 100
nSnps = 1000
pA1=0.5; pA2=0.5
epsilon=0.01
depths <- c(1:5,7.5,10,15)

for(D in c(0,0.05,0.1,0.25)){
  for(meanDepth in depths){

    pB1 <- 1-pA1
    pB2 <- 1-pA2

    ## Check that the contraint on D is met

    C1 <- round(max(-pA1*pA2,-pB1*pB2),16); C2 <- round(min(pB1*pA2,pA1*pB2),16)
    if (!((C1 <= D) & (D <= C2))) stop('Error: User specified value of D is outside it constrained region')

    ###### Calculate true values of the LD measures
    # r^2
    r2 <- D^2/(pA1*pA2*pB1*pB2)
    # D'
    Dcon <- ifelse(D<0, abs(C1), C2)
    Dprime <- D/(Dcon)

    ## Table of true haplotype probabilities
    haploP <- round(matrix(c(pA1*pA2+D,pA1*pB2-D,pB1*pA2-D,pB1*pB2+D), ncol=2,
                           dimnames=list(c('A1','B1'),c('A2','B2')), byrow=T),8)

    #Est.mat <- matrix(NA,nrow=length(Depths), ncol=10*3+2) # Matrix containing the monte carlo bias and variance of MLE's
    #colnames(Est.mat) <- c(paste0(rep(c("E(D","E(D'","E(r2","E(pa1","E(pa2","Se(D","Se(D'","Se(r2","Se(pa1","Se(pa2"),3),").",
    #                              c(rep(1,10),rep(2,10),rep(3,10))),"E(ep)","Se(ep)")

    trueGeno <- matrix(sapply(1:nSnps, function(x){
      obsH <- matrix(sample(c('AA','BA','AB','BB'),size=nInd*2,
                            prob=as.vector(haploP), replace=T), ncol=2)
      return(t(apply(obsH,1, function(x) c(sum(substr(x,1,1) == 'A'),sum(substr(x,2,2) == 'A')))))
    } ), nrow=nInd)
    depth <- matrix(rpois(nInd*2*nSnps,meanDepth),ncol=nSnps*2, nrow=nInd)

    ##### Generate the GBS data from the true genotypes
    # Simulate the observed A alleles
    aCounts <- matrix(rbinom(nInd*2*nSnps,depth,trueGeno/2),ncol=nSnps*2, nrow=nInd)
    bCounts <- depth - aCounts
    aCountsFinal <- matrix(rbinom(nInd*2*nSnps,aCounts,prob=1-epsilon),ncol=nSnps*2, nrow=nInd) + matrix(rbinom(nInd*2*nSnps,bCounts,prob=epsilon),ncol=nSnps*2, nrow=nInd)
    bCountsFinal <- depth - aCountsFinal
    #SEQgeno <- aCountsFinal/depth
    #SEQgeno[which(SEQgeno^2-SEQgeno<0)] <- 0.5
    #SEQgeno <- 2* SEQgeno  ## GBS genotype call
    #GBSfreq <- colMeans(SEQgeno,na.rm=TRUE)/2

    ## new way (which is like 14x faster)
    rd <- RA$new(list())
    temp <- UR$new(rd, ploid=2)
    temp$.__enclos_env__$private$updatePrivate(list(ref = aCountsFinal, alt=bCountsFinal, chrom=as.integer(1:(2*nSnps)), pos=as.integer(1:(2*nSnps)),
                                               nInd=as.integer(nInd), nSnps=as.integer(2*nSnps), gform="reference"))

    tt <- temp$.__enclos_env__$private$p_est(nClust = 3)
    temp$.__enclos_env__$private$updatePrivate(list(pfreq=tt[1,], ep=tt[2,]))
    GUSLD(temp,SNPpairs = cbind(seq(1,2*nSnps-1,2), seq(2,2*nSnps,2)), nClust=3, LDmeasure = c("Dcoef","Dprime","r2"), dp=12,
          file = paste0('test/GBS_quick_pA1_',pA1,'_pA2_',pA2,'_D',D,'_N',nInd,'_nSnps',nSnps,"_depth",meanDepth,'.txt',sep=''))
  }
}



## plot the results
D=c(0,0.05,0.15,0.25)

Est.mat <- Depths <- trueVar <- replicate(n = length(D), numeric(length(depths)), simplify=F)

for(i in 1:length(D)){
  for(j in 2:length(depths)){
    temp <- read.table(paste0('test/GBS_quick_pA1_',pA1,'_pA2_',pA2,'_D',D[i],'_N',nInd,'_nSnps',nSnps,"_depth",depths[j],'.txt_GUSLD.txt'),header=T,sep=",", stringsAsFactors = F)
    Est.mat[[i]][j] <- c(colMeans(temp[,1:3]))
    Depths[[i]] <- scan(paste0('SimData/GBS_depth_pA1_',pA1,'_pA2_',pA2,'_D',D[i],'_N',n,'Runs',Nruns,'.txt'))
    trueVar[[i]] <- scan(paste0('SimData/GBS_trueVar_pA1_',pA1,'_pA2_',pA2,'_D',D[i],'_N',n,'Runs',Nruns,'.txt'))
    ## remove depth 20
    Est.mat[[i]] <- Est.mat[[i]][-9,]
    Depths[[i]] <- Depths[[i]][-9]
  }
}







