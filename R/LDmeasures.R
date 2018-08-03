

Dcoef <- function(pA1,pA2,D){
  return(D)
}

Dprime <- function(pA1,pA2,D){
  C1 = max(-prod(c(pA1,pA2)),-prod(1-c(pA1,pA2))); C2 = min((1-pA1)*pA2,pA1*(1-pA2))
  return(D/ifelse(D<0,C1,C2))
}

r2 <- function(pA1,pA2,D){
  return(D^2/(prod(c(pA1,pA2,1-c(pA1,pA2)))))
}


