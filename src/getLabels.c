// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Boolean.h>
#include "utils.h"



SEXP getLabels(SEXP mod, SEXP dataset) {
	// compute responsibilities w.r.t classic EM style
	
	PROTECT(mod=coerceVector(mod, VECSXP));
	PROTECT(dataset=coerceVector(dataset, REALSXP));
	SEXP r_false = PROTECT(allocVector(INTSXP, 1));
	INTEGER(r_false)[0] = 0;	

	int k = length(coerceVector(getListElement(mod, "w"), REALSXP));
	int n = INTEGER(getAttrib(dataset, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(dataset, R_DimSymbol))[1];
	

	// build a  nxk matrix containing temp results
	gsl_matrix *response = gsl_matrix_calloc(n, k);
	gsl_vector_view view;
	
	SEXP current;
	for(int i=0; i<k; i++) {
		current = mvndensity(coerceVector(VECTOR_ELT(coerceVector(getListElement(mod, "mean"), VECSXP), i), REALSXP),
			coerceVector(VECTOR_ELT(coerceVector(getListElement(mod, "cov"), VECSXP), i), REALSXP), dataset, r_false);
		view = gsl_matrix_column(response, i);
		SXPtoVector(&(view.vector), current);
		gsl_vector_scale(&(view.vector), REAL(coerceVector(getListElement(mod, "w"), REALSXP))[i]);
	}
	
	gsl_vector *gslres = gsl_vector_alloc(n);
	
	for(int i=0; i<n; i++) {
		view = gsl_matrix_row(response, i);
		gsl_vector_set(gslres, i, (double)whichmax(&(view.vector))+1.0);
	}
	
	SEXP res;
	PROTECT(res=allocVector(INTSXP, n));
	intVectorToSXP(&res, gslres);
	
	gsl_matrix_free(response);
	gsl_vector_free(gslres);
	UNPROTECT(4);
	return(res);
	

}
