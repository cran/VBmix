#ifndef MKMEANS_H
#define MKMEANS_H

#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Boolean.h>
#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "utils.h"


SEXP mkmeans(SEXP, SEXP, SEXP, SEXP, SEXP);

#endif

