#ifndef EM_C
#define EM_C


#include "EM.h"

SEXP EM(SEXP data, SEXP ncomp, SEXP model, SEXP class, SEXP thres, SEXP maxit, SEXP rbic, 
	SEXP seeds, SEXP debug) {
	// FORTRAN and R-ext C code assume matrices in column major format.
	// M[i,j] is thus M[j*nr+i]

	GetRNGstate();

	// necessary if data is a matrix of int
	PROTECT(data = coerceVector(data, REALSXP));
	double *c_data = REAL(data);
	const char *c_model = CHAR(STRING_ELT(model,0));
	int n = INTEGER(getAttrib(data, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(data, R_DimSymbol))[1];
	int k = INTEGER(coerceVector(ncomp, INTSXP))[0];
	double c_thres = REAL(thres)[0];
	int c_maxit;
	if(maxit == R_NilValue) {
		c_maxit = INT_MAX;
	} else {
		c_maxit = INTEGER(coerceVector(maxit, INTSXP))[0];
	}
	int c_rbic = INTEGER(rbic)[0];
	int c_class = INTEGER(class)[0];
	int c_debug = INTEGER(debug)[0];
	int* c_seeds = INTEGER(seeds);
	
	if(strcmp(c_model, "general") && strcmp(c_model, "diagonal") && strcmp(c_model, "spherical")) {
		error("model should be either 'general', 'diagonal' or 'spherical'");
	}

	// to limit burden on the stack, everything that scales according to d (pot. large)
	// is allocated on the heap
	char *trans = "T";
	char *nontrans = "N";
	int first_tie = 2;
	double d_one =1.0;
	double alpha = 1.0;
	double d_zero=0.0;
	int i_one = 1;
	int i_zero = 0;
	int cells = n*d;
	int d2 = d*d;
	int dp1 = d+1;
	int info = 0;
	int ival;
	int cond1;
	int cond2;

	double dval;
	double dval2;
	double dvalk[k];

	double *eigvals;
	double *eigvecs;
	int *eigsupp;
	double *eigwork;
	int *eigiwork;
	if(!strcmp(c_model, "general")) {
		eigvals = calloc(d, sizeof(double));
		eigvecs = calloc(d*d, sizeof(double));
		eigsupp = calloc(2*d, sizeof(int));
		eigwork = calloc(26*d, sizeof(double));
		eigiwork = calloc(10*d, sizeof(int));
	}

	int eiglwork=26*d;
	int eigliwork=10*d;	

	double *c_mins = calloc(d, sizeof(double));
	double *c_maxs = calloc(d, sizeof(double));
	double *c_sampmean = calloc(d, sizeof(double));
	memset(c_mins, 0, d*sizeof(double));
	memset(c_maxs, 0, d*sizeof(double));
	for(int i=0; i<d; i++) {
		// use pointer arithmetics
		// warn: indexes as manipulated and returned by FORTRAN are 1-based !!!
		ival = imax(c_data+(i*n), n, 1);
		c_maxs[i] = c_data[i*n+ival];
		ival = imin(c_data+(i*n), n, 1);
		c_mins[i] = c_data[i*n+ival];
		c_sampmean[i] = sum(c_data+(i*n), n, 1) / n;
	}

	int noupdate[k];
	int nullcov[k];
	memset(noupdate, 0, k*sizeof(int));
	memset(nullcov, 0, k*sizeof(int));

	double c_w[k];
	double c_resp[k];
	double c_lnresp[k];
	double *c_samp = calloc(d, sizeof(double));
	int *c_labels = calloc(n, sizeof(int));
	double **c_mean = calloc(k, sizeof(double *));
	for(int i=0; i<k; i++) c_mean[i] = calloc(d, sizeof(double));
	double **c_cov = calloc(k, sizeof(double *));
	if(!strcmp(c_model, "general")) {
		for(int i=0; i<k; i++) c_cov[i] = calloc(d*d, sizeof(double));
	} else if(!strcmp(c_model, "diagonal")) {
		for(int i=0; i<k; i++) c_cov[i] = calloc(d, sizeof(double));
	} else {
		for(int i=0; i<k; i++) c_cov[i] = calloc(1, sizeof(double));
	}

	double *c_mincov;
	if(!strcmp(c_model, "general")) {
		c_mincov = calloc(d*d, sizeof(double));
	} else if(!strcmp(c_model, "diagonal")) {
		c_mincov = calloc(d, sizeof(double));
	} else {
		c_mincov = calloc(1, sizeof(double));
	}

	// temp accumulators
	double c_rnk[k];
	double **c_rnkx = calloc(k, sizeof(double *));
	for(int i=0; i<k; i++) c_rnkx[i] = calloc(d, sizeof(double));
	double **c_rnkxx = calloc(k, sizeof(double *));
	if(!strcmp(c_model, "general")) {
		for(int i=0; i<k; i++) c_rnkxx[i] = calloc(d*d, sizeof(double));
	} else if(!strcmp(c_model, "diagonal")) {
		for(int i=0; i<k; i++) c_rnkxx[i] = calloc(d, sizeof(double));
	} else {
		for(int i=0; i<k; i++) c_rnkxx[i] = calloc(1, sizeof(double));
	}

	double **c_inverse = calloc(k, sizeof(double *));
	if(!strcmp(c_model, "general")) {
		for(int i=0; i<k; i++) c_inverse[i] = calloc(d*d, sizeof(double));
	} else if(!strcmp(c_model, "diagonal")) {
		for(int i=0; i<k; i++) c_inverse[i] = calloc(d, sizeof(double));
	} else {
		for(int i=0; i<k; i++) c_inverse[i] = calloc(1, sizeof(double));
	}

	double c_lndetcov[k];

	for(int i=0; i<k; i++) c_w[i] = 1.0/k;
	memset(c_resp, 0, k*sizeof(double));
	memset(c_lnresp, 0, k*sizeof(double));
	memset(c_labels, 1, n*sizeof(int));
	for(int i=0; i<k; i++) {
		for(int j=0;j<d;j++) {
			//dval = c_mins[j] + (c_maxs[j] - c_mins[j]) * unif_rand();
			//c_mean[i][j] = dval;
			//c_mean[i][j] = c_sampmean[j] + (c_maxs[j] - c_mins[j]) * unif_rand() / 4.0;
			// k-means_like init
			c_mean[i][j] = c_data[j*n + (c_seeds[i]-1)];
			if(!strcmp(c_model, "general")) {
				c_cov[i][j*d+j] = pow(c_maxs[j] - c_mins[j], 2.0) / 4.0;
			} else if(!strcmp(c_model, "diagonal")) {
				c_cov[i][j] = pow(c_maxs[j] - c_mins[j], 2.0) / 4.0;
			} else {
				c_cov[i][0] += pow(c_maxs[j] - c_mins[j], 2.0) / 4.0;
			}
		}
		if(!strcmp(c_model, "spherical")) {
			c_cov[i][0] = c_cov[i][0] / d;
		}
	}
	
	for(int i=0; i<d; i++) {
		if(!strcmp(c_model, "general")) {
			c_mincov[i*d+i] = pow(c_maxs[i] - c_mins[i], 2.0) / 2048.0;
		} else if(!strcmp(c_model, "diagonal")) {
			c_mincov[i] = pow(c_maxs[i] - c_mins[i], 2.0) / 2048.0;
		} else {
			c_mincov[0] += pow(c_maxs[i] - c_mins[i], 2.0) / 2048.0;
		}
	}
	if(!strcmp(c_model, "spherical")) {
		c_mincov[0] = c_mincov[0] / d;
	}

	int done=0;
	int it = 1;
	double oldLikelihood = -DBL_MAX;
	double newLikelihood;
	double rLikelihood;
	double c_bic;

	// init covariance decompositions
	for(int i=0;i<k;i++) symdecomp(c_cov[i], c_inverse[i], c_lndetcov+i, d, c_model);


	while(!done) {

		// reset sufficient statistics
		memset(c_rnk, 0, k*sizeof(double));
		for(int i=0;i<k;i++) memset(c_rnkx[i], 0, d*sizeof(double));
		if(!strcmp(c_model, "general")) {
			for(int i=0;i<k;i++) memset(c_rnkxx[i], 0, d*d*sizeof(double));
		} else if(!strcmp(c_model, "diagonal")) {
			for(int i=0;i<k;i++) memset(c_rnkxx[i], 0, d*sizeof(double));
		} else {
			for(int i=0;i<k;i++) memset(c_rnkxx[i], 0, 1*sizeof(double));
		}

		// loop over elements
		for(int i=0; i<n; i++) {
			// compute resp (E step) and update likelihood inline
			F77_CALL(dcopy)(&d, c_data+i, &n, c_samp, &i_one);
			for(int j=0; j<k; j++) {
				
				c_resp[j] = c_w[j] * dmnorm(c_samp, c_mean[j], c_lndetcov[j], c_inverse[j], d, 0, c_model);
				if(c_w[j] > 0) {
					c_lnresp[j] = log(c_w[j]) + dmnorm(c_samp, c_mean[j], c_lndetcov[j], c_inverse[j], d, 1, c_model);
				} else {
					c_lnresp[j] = -DBL_MAX;
				}
			}

			// compute new labels
			dval = sum(c_resp, k, 1);
			if(dval > 0 && c_class != TRUE) {
				dval = 1.0/dval;
				F77_CALL(dscal)(&k, &dval, c_resp, &i_one);
			} else {
				memset(c_resp, 0, k*sizeof(double));
				c_resp[imax(c_lnresp, k, 1)] = 1.0;
			}
			c_labels[i] = imax(c_resp, k, 1)+1;
			

			// compute sufficient statistics
			alpha = 1.0;
			F77_CALL(daxpy)(&k, &alpha, c_resp, &i_one, c_rnk, &i_one);
			for(int j=0; j<k; j++) {
				F77_CALL(dcopy)(&d, c_data+i, &n, c_samp, &i_one);
				F77_CALL(daxpy)(&d, c_resp+j, c_samp, &i_one, c_rnkx[j], &i_one);
				if(!strcmp(c_model, "general")) {
					F77_CALL(dsyr)("L", &d, c_resp+j, c_samp, &i_one, c_rnkxx[j], &d);
				} else if(!strcmp(c_model, "diagonal")) {
					for(int m=0; m<d; m++) c_rnkxx[j][m] += c_resp[j] * pow(c_samp[m], 2.0);
				} else {
					for(int m=0; m<d; m++) c_rnkxx[j][0] += c_resp[j] * pow(c_samp[m], 2.0);
				}
			}

		}

		if(c_debug == TRUE) {
			for(int i=0; i<k; i++) Rprintf("%5.1f ", c_rnk[i]);
			Rprintf("\n");
		}

		// update model (M step)
		for(int i=0; i<k; i++) {
			if(c_rnk[i] < 2.0 && !noupdate[i]) {
				// zero components with no sufficient support
				c_rnk[i] = 0.0;
				noupdate[i] = 1;
				oldLikelihood = -DBL_MAX;
			} else if(!noupdate[i]) {
				dval = 1.0/(c_rnk[i]);
				F77_CALL(dscal)(&d, &dval, c_rnkx[i], &i_one);
				if(!strcmp(c_model, "general")) {
					F77_CALL(dscal)(&d2, &dval, c_rnkxx[i], &i_one);
				} else if(!strcmp(c_model, "diagonal")) {
					F77_CALL(dscal)(&d, &dval, c_rnkxx[i], &i_one);
				} else {
					c_rnkxx[i][0] = c_rnkxx[i][0] * dval;
				}
				
				F77_CALL(dcopy)(&d, c_rnkx[i], &i_one, c_mean[i], &i_one);

				// do not process covariance of components that previously diverged
				// towards singularities
				if(!nullcov[i]) {
					if(!strcmp(c_model, "general")) {
						F77_CALL(dcopy)(&d2, c_rnkxx[i], &i_one, c_cov[i], &i_one);
						alpha = -1.0;
						F77_CALL(dsyr)("L", &d, &alpha, c_mean[i], &i_one, c_cov[i], &d);
					} else if(!strcmp(c_model, "diagonal")) {
						for(int j=0; j<d; j++) c_cov[i][j] = c_rnkxx[i][j] - pow(c_mean[i][j], 2.0);
					} else {
						c_cov[i][0] = c_rnkxx[i][0];
						for(int j=0; j<d; j++) c_cov[i][0] -= pow(c_mean[i][j], 2.0);
						c_cov[i][0] = c_cov[i][0] / d;
					}
			

					// test singularity of components (mostly based on eigenvals)
					if(!strcmp(c_model, "general")) {
						alpha = 0.0;
						// dsyevr modifies diagonalized matrix, so use a proxy
						F77_CALL(dcopy)(&d2, c_cov[i], &i_one, c_rnkxx[i], &i_one);
						F77_CALL(dsyevr)("N", "A", "L", &d, c_rnkxx[i], &d, &alpha, &alpha, &i_one, &i_one,
							&alpha, &ival, eigvals, eigvecs, &d, eigsupp, eigwork, &eiglwork, eigiwork,
							&eigliwork, &info);
						//dval = 0.0;
						//for(int j=0; j<d; j++) dval += c_mincov[j*d+j];
						dval = pow(10.0, -6.0);
						cond1 = (sum(eigvals, d, 1) < dval);
						// eigenvalues in ascending order
						// use them to rescale component gracefully
						cond2 = (eigvals[0] / eigvals[d-1] < pow(10.0, -6.0));
						if(cond1 || cond2) {
							if(c_debug == TRUE) Rprintf("comp %d cond1 %d cond2 %d\n", i, cond1, cond2);
							nullcov[i] = 1;
							oldLikelihood = -DBL_MAX;
							memset(c_cov[i], 0, d*d*sizeof(double));
							for(int j=0; j<d; j++) c_cov[i][j*d+j] = max(dval, sum(eigvals, d, 1));
						}

					} else if(!strcmp(c_model, "diagonal")) {
						// same block adapted to diagonal matrix
						ival = imin(c_cov[i], d, 1);
						dval = c_cov[i][ival];
						ival = imax(c_cov[i], d, 1);
						dval = dval / c_cov[i][ival];
						dval2 = sum(c_mincov, d, 1);
						cond1 = (sum(c_cov[i], d, 1) < dval2);
						cond2 = (dval < pow(10.0, -6.0));
						if(cond1 || cond2) {
							nullcov[i] = 1;
							oldLikelihood = -DBL_MAX;
							dval = sum(c_cov[i], d, 1);
							for(int j=0; j<d; j++) c_cov[i][j] = max(dval, dval2);
						}
					} else {
						if(c_cov[i][0] < c_mincov[0]) {
							c_cov[i][0] = c_mincov[0];
							nullcov[i] = 1;
							oldLikelihood = -DBL_MAX;								
						}
					}
				}
	
			}
		}

		// check sums/weights vector, and normalize
		dval = sum(c_rnk, k, 1);
		// in theory dval should be n, but this moight not be true in case of singularities
		if(dval == 0.0) {
			for(int j=0;j<k; j++) c_rnk[j] = 1.0/k;
			dval = 1.0;
		}
		dval = 1.0/dval;
		F77_CALL(dcopy)(&k, c_rnk, &i_one, c_w, &i_one);
		F77_CALL(dscal)(&k, &dval, c_w, &i_one);

		// restore symmetric cov matrices
		if(!strcmp(c_model, "general")) {
			for(int i=0;i<k;i++) lowerComplete(c_cov[i], d);
		}
		

		// update cov decompositions
		for(int i=0;i<k;i++) symdecomp(c_cov[i], c_inverse[i], c_lndetcov+i, d, c_model);

		// update likelihood and bic, test on continuation
		newLikelihood = 0.0;
		rLikelihood = 0.0;
		for(int i=0; i<n; i++) {
			F77_CALL(dcopy)(&d, c_data+i, &n, c_samp, &i_one);
			dval=0.0;
			for(int j=0; j<k; j++) {
				dval += c_w[j] * dmnorm(c_samp, c_mean[j], c_lndetcov[j], c_inverse[j], d, 0, c_model);
			}
			newLikelihood += log(dval);
			dval=0.0;
			for(int j=0; j<k; j++) {
				dval += c_w[j] * drmnorm(c_samp, c_mean[j], c_lndetcov[j], c_inverse[j], d, 0, c_model);
			}
			rLikelihood += log(dval);
		}


		it += 1;
		dval = newLikelihood - oldLikelihood;

		if(c_debug == TRUE) {
			char str1[9];
			char str2[9];
			if(oldLikelihood == -DBL_MAX) {
				sprintf(str1, "-Inf");
			} else {
				sprintf(str1, "%9.2f", oldLikelihood);
			} 
			if(newLikelihood == -DBL_MAX) {
				sprintf(str2, "-Inf");
			} else {
				sprintf(str2, "%9.2f", newLikelihood);
			}
			Rprintf("old: %s, new: %s ", str1, str2);
		}

		if(dval < c_thres || it > c_maxit) {
			done = 1;
			if(c_debug == TRUE) Rprintf("done!");
		}

		if(c_debug == TRUE) Rprintf("\n");

		oldLikelihood = newLikelihood;
		if(c_rbic == TRUE) {
			c_bic = -2.0 * rLikelihood;
		} else {
			c_bic = -2.0 * newLikelihood;
		}
		
		if(!strcmp(c_model, "spherical")) {
			c_bic += (k-1 + k*(d+1)) * log(n);
		} else if(!strcmp(c_model, "diagonal")) {
			c_bic += (k-1 + 2*k*d) * log(n);
		} else {
			c_bic += (k-1 + k*d + k*(d*(d+1))/2.0) * log(n);
		}


	}






	// formatting R output
	SEXP res, w, mean, cov, labels, likelihood, bic, cur, names;

	PROTECT(res = allocVector(VECSXP, 6));
	PROTECT(w = allocVector(REALSXP, k));
	PROTECT(mean = allocVector(VECSXP, k));
	PROTECT(cov = allocVector(VECSXP, k));
	PROTECT(labels = allocVector(INTSXP, n));
	PROTECT(likelihood = allocVector(REALSXP, 1));
	PROTECT(bic = allocVector(REALSXP, 1));

	F77_CALL(dcopy)(&k, c_w, &i_one, REAL(w), &i_one);
	for(int i=0; i<n; i++) INTEGER(labels)[i] = c_labels[i];
	if(c_rbic == TRUE) {
		REAL(likelihood)[0] = rLikelihood;
	} else {
		REAL(likelihood)[0] = newLikelihood;
	}
	REAL(bic)[0] = c_bic;


	for(int i=0; i<k; i++) {
		PROTECT(cur=allocVector(REALSXP, d));
		F77_CALL(dcopy)(&d, c_mean[i], &i_one, REAL(cur), &i_one);
		SET_VECTOR_ELT(mean, i, cur);
		UNPROTECT(1);
		// restore full matrices on output
		PROTECT(cur=allocMatrix(REALSXP, d, d));
		memset(REAL(cur), 0, d*d*sizeof(double));
		if(!strcmp(c_model, "general")) {
			F77_CALL(dcopy)(&d2, c_cov[i], &i_one, REAL(cur), &i_one);
		} else if(!strcmp(c_model, "diagonal")) {
			F77_CALL(dcopy)(&d, c_cov[i], &i_one, REAL(cur), &dp1);
		} else {
			for(int j=0; j<d; j++) REAL(cur)[j*d+j] = c_cov[i][0];
		}

		SET_VECTOR_ELT(cov, i, cur);
		UNPROTECT(1);	
	}
	
	if(!strcmp(c_model, "general")) {
		free(eigvals);
		free(eigvecs);
		free(eigsupp);
		free(eigwork);
		free(eigiwork);
	}
	free(c_mins);
	free(c_maxs);
	free(c_samp);
	free(c_labels);
	free(c_mincov);
	for(int i=0;i<k;i++) {
		free(c_mean[i]);
		free(c_cov[i]);
		free(c_rnkx[i]);
		free(c_rnkxx[i]);
		free(c_inverse[i]);
	}
	free(c_mean);
	free(c_cov);
	free(c_rnkx);
	free(c_rnkxx);
	free(c_inverse);


	SET_VECTOR_ELT(res, 0, w);
	SET_VECTOR_ELT(res, 1, mean);
	SET_VECTOR_ELT(res, 2, cov);
	SET_VECTOR_ELT(res, 3, labels);
	SET_VECTOR_ELT(res, 4, likelihood);
	SET_VECTOR_ELT(res, 5, bic);


	PROTECT(names = allocVector(VECSXP, 6));
	SET_VECTOR_ELT(names, 0, mkChar("w"));
	SET_VECTOR_ELT(names, 1, mkChar("mean"));
	SET_VECTOR_ELT(names, 2, mkChar("cov"));
	SET_VECTOR_ELT(names, 3, mkChar("labels"));
	SET_VECTOR_ELT(names, 4, mkChar("likelihood"));
	SET_VECTOR_ELT(names, 5, mkChar("bic"));
	setAttrib(res, R_NamesSymbol, names);

	UNPROTECT(9);
	PutRNGstate();
	return(res);
	//return(R_NilValue);
}

#endif

