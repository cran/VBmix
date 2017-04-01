#ifndef MKMEANS_C
#define MKMEANS_C


#include "mkmeans.h"

// as at each step a constraint is set on covariance matrices,
// and responsibilities are discretized, should bear some links with CEM - C. Schmidt
// approaches. The constraint can be viewed as a way to handle the tendency of EM 
// to "absorb" all the data distribution in one single component (with single example)
// when here we want to reach good classification models on potentially inorthodox
// distributions
SEXP mkmeans(SEXP dat, SEXP k, SEXP maxiter, SEXP seeds, SEXP prior) {

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
	SEXP r_true = PROTECT(allocVector(INTSXP, 1));
	SEXP r_false = PROTECT(allocVector(INTSXP, 1));
	INTEGER(r_true)[0] = 1;
	INTEGER(r_false)[0] = 0;

	int* c_seeds = INTEGER(coerceVector(seeds, INTSXP));
	int c_prior = INTEGER(coerceVector(prior, INTSXP))[0];

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


	
	for(int i=0; i<d; i++) {
		c_mins[i] = vmin(c_dat+(i*n), n, 1);
		c_maxs[i] = vmax(c_dat+(i*n), n, 1);
		d_coeff = max(pow(c_maxs[i] - c_mins[i], 2.0), 8*DBL_EPSILON);
		c_regcov[i] = d_coeff / 4.0;

		//c_regcov[i] = (c_maxs[i] - c_mins[i]) / 16.0;
		// penalizes smaller clusters
	}


	double *c_centers = REAL(centers);
	double **c_covs = calloc(c_k, sizeof(double *));
	for(int i=0; i<c_k; i++) {
		F77_CALL(dcopy)(&d, c_dat+(c_seeds[i]), &n, c_centers+i, &c_k);
		c_covs[i] = REAL(VECTOR_ELT(covs, i));
		memset(c_covs[i], 0, d*d*sizeof(double));
		for(int j=0; j<d; j++) {
			//c_covs[i][j*(d+1)] = 1.0;
			c_covs[i][j*(d+1)] = c_regcov[j];
		}
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

	// for matrix rescale
	double avg_extent=0.0;
	double extents[c_k];
	double evals[c_k][d];
	double evecs[c_k][d2];
	double w_tridiag[d-1];
	double w_tau[d-1];
	int w_lwork = d*3;
	double w_work[w_lwork];
	int w_info = 0;



	for(nit=0; nit<c_maxit; nit++) {
		changed = FALSE;
		for(int i=0; i<c_k; i++) sumlabs[i] = 0;
		dists = PROTECT(gdist(dat, centers, covs, r_false));
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
				//if(nit % 5 == 0) {
					memset(c_covs[i], 0, d*d*sizeof(double));
				//}
				if(sumlabs[i] > 0) {
					for(int j=0; j<d; j++) c_centers[i+j*c_k] = 0.0;
				} 
				// avoid biasing average step
				//if(sumlabs[i] < 2) {
				//	for(int j=0; j<d; j++) c_covs[i][j*(d+1)] = 1.0;
				//}
			}

			for(int i=0; i<n; i++) {
				d_coeff = 1.0;
				F77_CALL(daxpy)(&d, &d_coeff, c_dat+i, &n, c_centers+c_newlabs[i], &c_k);
			}

			for(int i=0; i<c_k; i++) {
				if(sumlabs[i] > 0) {
					d_coeff = 1.0/(double)(sumlabs[i]);
					F77_CALL(dscal)(&d, &d_coeff, c_centers+i, &c_k);
				}
			}

			//if(nit % 5 == 0) {
				for(int i=0; i<n; i++) {
					F77_CALL(dcopy)(&d, c_dat+i, &n, d_vec, &i_one);
					d_coeff = -1.0;
					F77_CALL(daxpy)(&d, &d_coeff, c_centers+c_newlabs[i], &c_k, d_vec, &i_one);
					d_coeff = 1.0;
					if(sumlabs[c_newlabs[i]] > 1) {
						F77_CALL(dsyr)("L", &d, &d_coeff, d_vec, &i_one, c_covs[c_newlabs[i]], &d);
					}
				}				
			//}



//#     # rescale so that geom. mean of eigenvals is res.extent.ext[i]
//#     eig <- eigen(Sigma)
//#     V <- eig$vectors
//#     L <- eig$values
//#     L <- L * res.extent.ext[i] / geom.mean(L)
//#     Sigma <- V %*% diag(L) %*% t(V)


			double avg_trace=0.0;
			double traces[c_k];
			memset(traces, 0, c_k*sizeof(double));
			for(int i=0; i<c_k; i++) {
				//if(nit % 5 == 0) {
					for(int j=0; j<d; j++) {
						c_covs[i][j*(d+1)] += c_prior*c_regcov[j];
					}

					if(sumlabs[i] > 1) {
						d_coeff = 1.0/((double)(sumlabs[i]+c_prior-1));
						F77_CALL(dscal)(&d2, &d_coeff, c_covs[i], &i_one);
					}
				//}
				for(int j=0; j<d; j++) traces[i] += c_covs[i][j*(d+1)];
				avg_trace += traces[i];
			}

			avg_trace = avg_trace / (double)c_k;
			for(int i=0; i<c_k; i++) {
				// to encourage balanced classes, rescale so that trace is d
				// PBR : not scale independent then... better rescale so that trace is
				// avg-traces, ie avg trace remains the same.
				// somehow enforcement of equality constraints
				d_coeff = avg_trace / traces[i];
				F77_CALL(dscal)(&d2, &d_coeff, c_covs[i], &i_one);
			}



//			avg_extent=0.0;
//			w_info = 0;
//
////			for(int i=0; i<c_k; i++) {
////				memset(evals[i], 0, d*sizeof(double));
////				memset(evecs[i], 0, d2*sizeof(double));
////			}
//
//			for(int i=0; i<c_k; i++) {
//				if(sumlabs[i] > 1) {
//					d_coeff = 1.0/(double)(sumlabs[i]-1);
//					F77_CALL(dscal)(&d2, &d_coeff, c_covs[i], &i_one);
//				}
//				for(int j=0; j<d; j++) {
//					c_covs[i][j*(d+1)] += c_regcov[j];
//				}
//
//				//F77_CALL(dsytrd)("L", &d, c_covs[i], &d, evals[i], w_tridiag, &w_tau, w_work, &w_lwork, &w_info);
//				//F77_CALL(dsteqr)("V", &d, &w_diag, &w_tridiag, evecs[i], &d, &w_lwork, &w_info);
//				F77_CALL(dcopy)(&d2, c_covs[i], &i_one, evecs[i], &i_one);
//				F77_CALL(dsyev)("V", "L", &d, evecs[i], &d, evals[i], w_work, &w_lwork, &w_info);
//				extents[i]=1.0;
//				for(int j=0; j<d; j++) {
//					extents[i] *= evals[i][j];
//				}
//				extents[i] = pow(extents[i], 1.0/(double)d);
//				avg_extent += extents[i];
//			}
//
//			avg_extent = avg_extent / (double)c_k;
//
//			// according to matrix properties, U \Lambda U^T boild down to
//			// \sum \lambda_i u_i u_i^T
//			for(int i=0; i<c_k; i++) {
//				d_coeff = avg_extent / extents[i];
//				F77_CALL(dscal)(&d, &d_coeff, evals[i], &i_one);
//				memset(c_covs[i], 0, d2*sizeof(double));
//				for(int j=0; j<d; j++) {
//					F77_CALL(dsyr)("L", &d, &(evals[i][j]), evecs[i]+(i*d), &i_one, c_covs[i], &d);
//				}
//			}


//			for(int i=0; i<c_k; i++) {
//				for(int j=0; j<d; j++) {
//					c_covs[i][j*(d+1)] += c_regcov[j];
//				}
//				// decayingly less importance of the regularization
//				if(sumlabs[i] > 1) {
//					d_coeff = 1.0/(double)(sumlabs[i]-1);
//					F77_CALL(dscal)(&d2, &d_coeff, c_covs[i], &i_one);
//				}
//			}


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

	UNPROTECT(10);
	//UNPROTECT(8);
	PutRNGstate();
	return(res);
}

#endif

