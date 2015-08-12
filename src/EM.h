#ifndef EM_H
#define EM_H

#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Boolean.h>
#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "utils.h"


SEXP EM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif

