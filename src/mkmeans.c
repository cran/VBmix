#ifndef MKMEANS_C
#define MKMEANS_C


#include "mkmeans.h"

SEXP mkmeans(SEXP dat, SEXP k, SEXP maxiter, SEXP seeds) {

	GetRNGstate();

	// necessary if data is a matrix of int
	dat = PROTECT(coerceVector(dat, REALSXP));
	double *c_dat = REAL(dat);
	int n = INTEGER(getAttrib(dat, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(dat, R_DimSymbol))[1];
	int c_k = INTEGER(coerceVector(k, INTSXP))[0];
	int c_maxit = INTEGER(coerceVector(maxiter, INTSXP))[0];
	int i_one = 1;
	double d_coeff;
	int d2 = d*d;
	double *d_vec = calloc(d, sizeof(double));

	int* c_seeds = INTEGER(seeds);

	double *c_mins = calloc(d, sizeof(double));
	double *c_maxs = calloc(d, sizeof(double));
	double *c_regcov = calloc(d, sizeof(double));
	memset(c_mins, 0, d*sizeof(double));
	memset(c_maxs, 0, d*sizeof(double));
	memset(c_regcov, 0, d*sizeof(double));

	// params to gdist established as SEXP data
	SEXP centers = PROTECT(allocMatrix(REALSXP, c_k, d));
	SEXP covs = PROTECT(allocVector(VECSXP, c_k));
	SEXP cur;
	for(int i=0; i<c_k; i++) {
		SET_VECTOR_ELT(covs, i, allocMatrix(REALSXP, d, d));
	}

	double *c_centers = REAL(centers);
	double **c_covs = calloc(c_k, sizeof(double *));
	for(int i=0; i<c_k; i++) {
		F77_CALL(dcopy)(&d, c_dat+(c_seeds[i]), &n, c_centers+i, &c_k);
		c_covs[i] = REAL(VECTOR_ELT(covs, i));
		memset(c_covs[i], 0, d*d*sizeof(double));
		for(int j=0; j<d; j++) {
			c_covs[i][j*(d+1)] = 1.0;
		}
	}
	
	for(int i=0; i<d; i++) {
		c_mins[i] = vmin(c_dat+(i*n), n, 1);
		c_maxs[i] = vmax(c_dat+(i*n), n, 1);
		c_regcov[i] = pow(c_maxs[i] - c_mins[i], 2.0) / 2048.0;
	}


	SEXP newlabs = PROTECT(allocVector(INTSXP, n));
	int *c_newlabs = INTEGER(newlabs);
	for(int i=0; i<n; i++) {
		c_newlabs[i] = 0;
	}
	SEXP dists;
	double *c_dists;
	int nit;
	int sumlabs[c_k];
	Rboolean changed;

	for(nit=0; nit<c_maxit; nit++) {
		changed = FALSE;
		for(int i=0; i<c_k; i++) sumlabs[i] = 0;
		dists = PROTECT(gdist(dat, centers, covs));
		c_dists = REAL(dists);

		for(int i=0; i<n; i++) {
			int newlab = imin(c_dists+i, c_k, n);
			sumlabs[newlab]++;
			if(c_newlabs[i] != newlab) {
				c_newlabs[i] = newlab;
				changed = TRUE;
			}
		}

		UNPROTECT(1);

		if(changed) {
			for(int i=0; i<c_k; i++) {
				memset(c_covs[i], 0, d*d*sizeof(double));
				if(sumlabs[i] > 0) {
					for(int j=0; j<d; j++) c_centers[i+j*c_k] = 0.0;
				} 
				if(sumlabs[i] < 2) {
					for(int j=0; j<d; j++) c_covs[i][j*(d+1)] = 1.0;
				}
			}

			for(int i=0; i<n; i++) {
				d_coeff = 1.0;
				F77_CALL(daxpy)(&d, &d_coeff, c_dat+i, &n, c_centers+c_newlabs[i], &c_k);
			}

			for(int i=0; i<c_k; i++) {
				d_coeff = 1.0/(double)(sumlabs[i]);
				F77_CALL(dscal)(&d, &d_coeff, c_centers+i, &c_k);
			}

			for(int i=0; i<n; i++) {
				F77_CALL(dcopy)(&d, c_dat+i, &n, d_vec, &i_one);
				d_coeff = -1.0;
				F77_CALL(daxpy)(&d, &d_coeff, c_centers+c_newlabs[i], &c_k, d_vec, &i_one);
				d_coeff = 1.0;
				if(sumlabs[c_newlabs[i]] > 1) {
					F77_CALL(dsyr)("L", &d, &d_coeff, d_vec, &i_one, c_covs[c_newlabs[i]], &d);
				}
			}

			for(int i=0; i<c_k; i++) {
				if(sumlabs[i] > 1) {
					d_coeff = 1.0/(double)(sumlabs[i]-1);
					F77_CALL(dscal)(&d2, &d_coeff, c_covs[i], &i_one);
				}

				d_coeff = 0.0;
				for(int j=0; j<d; j++) {
					c_covs[i][j*(d+1)] += c_regcov[j];
					d_coeff += c_covs[i][j*(d+1)];
				}
				// to encourage balanced classes, rescale so that trace is d 
				d_coeff = (double)d / d_coeff;
				F77_CALL(dscal)(&d2, &d_coeff, c_covs[i], &i_one);				
			}


		}
	}


	// only eventually
	for(int i=0; i<c_k; i++) fillUpper(c_covs[i], d);

	// prepare output
	SEXP res = PROTECT(allocVector(VECSXP, 4));
	SEXP w = PROTECT(allocVector(REALSXP, c_k));
	SEXP mean = PROTECT(allocVector(VECSXP, c_k));

	for(int i=0; i<c_k; i++) {
		REAL(w)[i] = (double)(sumlabs[i]) / (double)n;
		SET_VECTOR_ELT(mean, i, allocVector(REALSXP, d));
		F77_CALL(dcopy)(&d, c_centers+i, &c_k, REAL(VECTOR_ELT(mean, i)), &i_one);
	}

	// add 1 to all labels for output
	for(int i=0; i<n; i++) c_newlabs[i] += 1;

	SET_VECTOR_ELT(res, 0, newlabs);
	SET_VECTOR_ELT(res, 1, w);
	SET_VECTOR_ELT(res, 2, mean);
	SET_VECTOR_ELT(res, 3, covs);

	SEXP names = PROTECT(allocVector(VECSXP, 4));
	SET_VECTOR_ELT(names, 0, mkChar("labels"));
	SET_VECTOR_ELT(names, 1, mkChar("w"));
	SET_VECTOR_ELT(names, 2, mkChar("mean"));
	SET_VECTOR_ELT(names, 3, mkChar("cov"));
	setAttrib(res, R_NamesSymbol, names);

	free(d_vec);
	free(c_mins);
	free(c_maxs);
	free(c_regcov);
	free(c_covs);

	//UNPROTECT(9);
	UNPROTECT(8);
	PutRNGstate();
	return(res);
}

#endif

