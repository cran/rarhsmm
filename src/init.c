//
//  init.c
//  
//
//  Created by XuZekun on 3/3/17.
//
//

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP rarhsmm_rmultinomial(SEXP,SEXP,SEXP,SEXP);
extern SEXP rarhsmm_mvrnorm(SEXP,SEXP,SEXP);
extern SEXP rarhsmm_mvdnorm(SEXP,SEXP,SEXP,SEXP);

//extern SEXP RcppArmadillo_armadillo_version(SEXP);
//extern SEXP RcppArmadillo_fastLm(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rarhsmm_rmultinomial", (DL_FUNC) &rarhsmm_rmultinomial, 4}, //number of parms
    {"rarhsmm_mvrnorm", (DL_FUNC) &rarhsmm_mvrnorm, 3},
    {"rarhsmm_mvdnorm", (DL_FUNC) &rarhsmm_mvdnorm, 4},
    {NULL, NULL, 0}
};

void R_init_rarhsmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
