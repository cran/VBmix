// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <R.h>
#include <Rinternals.h>

SEXP couple(SEXP tab1, SEXP tab2) {
	
	SEXP res;
	double val=0.0;
	int n, i, j;
	
	PROTECT(tab1 = coerceVector(tab1, INTSXP));
	PROTECT(tab2 = coerceVector(tab2, INTSXP));
	
	int *ptr1 = INTEGER(tab1);
	int *ptr2 = INTEGER(tab2);
	double *ptr;
	
	// assume tab1 and tab2 are same length
	n = length(tab1);
	
	for(i=0; i<(n-1); i++) {
		for(j=(i+1); j<n; j++) {
			if(ptr1[i] == ptr1[j]) {
				if(ptr2[i] != ptr2[j]) {
					val = val + 1.0;
				}
			} else {
				if(ptr2[i] == ptr2[j]) {
					val = val + 1.0;
				}
			}
		}
	}
	
	val = val * 2.0 / (n * (n-1));
	
	PROTECT(res = allocVector(REALSXP, 1));
	ptr = REAL(res);
	ptr[0] = val;
	
	UNPROTECT(3);
	return(res);
}
