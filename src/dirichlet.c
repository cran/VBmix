// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_randist.h>

SEXP dDirichlet(SEXP s_alpha, SEXP s_x1, SEXP s_x2) {
	// return the dirichlet lnpdf
	double alpha = REAL(s_alpha)[0];
	double x1 = REAL(s_x1)[0];
	double x2 = REAL(s_x2)[0];
	
	SEXP result;
	PROTECT(result=allocVector(REALSXP, 1));
	
	double alphatab[3];
	double xtab[3];
	
	for(int i=0; i<3; i++) alphatab[i]=alpha;
	xtab[0] = x1;
	xtab[1] = x2;
	xtab[2] = 1.0-x1-x2;
	
	REAL(result)[0] = gsl_ran_dirichlet_pdf(3, alphatab, xtab);
	UNPROTECT(1);
	return(result);
}
