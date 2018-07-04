
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP ll_gusld_c(SEXP LD, SEXP p, SEXP ep, SEXP AA, SEXP AB, SEXP BB, SEXP nInd){

  int ind, ind2, nInd_c, *pAA, *pAB, *pBB;
  double sum, LD_c;
  nInd_c = INTEGER(nInd)[0];
  double ep_c[2] = {REAL(ep)[0], REAL(ep)[1]} ;
  double p_c[2] = {REAL(p)[0], REAL(p)[0]};
  LD_c = REAL(LD)[0];
  // Define the pointers to the other input R variables
  pAA = INTEGER(AA);
  pAB = INTEGER(AB);
  pAB = INTEGER(BB);

  // compute the haplotype probabilities for given parameter values
  double h1 = p_c[0]*p_c[1] + LD_c;
  double h2 = p_c[0]*(1-p_c[1]) - LD_c;
  double h3 = (1-p_c[0])*p_c[1] - LD_c;
  double h4 = (1-p_c[0])*(1-p_c[1]) + LD_c;
  double p11 = h1*h1;
  double p21 = 2*h1*h3;
  double p31 = h3*h3;
  double p12 = 2*h1*h2;
  double p32 = 2*h3*h4;
  double p13 = h2*h2;
  double p23 = 2*h2*h4;
  double p33 = h4*h4;
  double p22 = 2*h1*h4+2*h2*h3;
  // Set up the variables to be returned
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  double llval = 0, *pll;
  pll = REAL(ll);

  // Compute the likelihood value
  for(ind = 0; ind < nInd_c; ind++){
    sum = 0;
    ind2 = nInd_c+ind;
    llval += log(pAA[ind]*(pAA[ind2]*p11 + pAB[ind2]*p12 + pBB[ind2]*p13) +
      pAB[ind]*(pAA[ind2]*p21 + pAB[ind2]*p22 + pBB[ind2]*p23) +
      pBB[ind]*(pAA[ind2]*p31 + pAB[ind2]*p32 + pBB[ind2]*p33));
  }

  // return the likelihood value
  pll[0] = llval;
  UNPROTECT(1);
  return ll;
}



