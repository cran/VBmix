// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "utils.h"


static double max_val(double *vec, int size) {
	int i;
	double themax = GSL_NEGINF;
	for(i=0; i<size; i++) {
		themax = (vec[i] > themax) ? vec[i] : themax;
	}
	return(themax);
}

	

SEXP vbcomp(SEXP models, SEXP ncomp, SEXP thres, SEXP maxit) {
	// for convenience get pointers to data structures
	SEXP w, mean, cov;
	w = getListElement(models, "w");
	mean = getListElement(models, "mean");
	cov = getListElement(models, "cov");
	// init prior model + precalc quantities (inverse wishart, constant terms in bound)
	// get sizes from data
	int L = length(w);
	int k = INTEGER(coerceVector(ncomp, INTSXP))[0];
	int d = length(VECTOR_ELT(mean, 0));
	int N = 200000;
	
	
	// declare stack variables for all usages
	int i, j, u, v, ind;
	double interm1;
	double interm2;

	double *read;

	gsl_vector_view view1;
	gsl_vector_view view2;
	gsl_vector_view view3;
	gsl_matrix_view matview1;
	gsl_matrix_view matview2;
	gsl_matrix_view matview3;
	double tempmat1[d*d];
	double tempmat2[d*d];
	double tempvec1[d];
	double tempvec2[d];
	gsl_permutation *perm = gsl_permutation_alloc(d);
	

	
	
	// first get maxs and mins from the models
	double mins[d];
	double maxs[d];
	
	for(i=0; i<d; i++) {
		mins[i] = GSL_POSINF;
		maxs[i] = GSL_NEGINF;
	}
	
	for(i=0; i<L; i++) {
		view1 = gsl_vector_view_array(REAL(VECTOR_ELT(mean, i)), d);
		for(j=0; j<d; j++) {
			if(gsl_vector_get(&view1.vector, j) > maxs[j]) {
				maxs[j] = gsl_vector_get(&view1.vector, j);
			}
			if(gsl_vector_get(&view1.vector, j) < mins[j]) {
				mins[j] = gsl_vector_get(&view1.vector, j);
			}
		}
	}
	
	// init estim model with prior model (+ ln det wishart + inv wish_0 + ln C alpha)
	// only prior means with k values : all other are identical for all k.
	double alpha_prior;
	double beta_prior;
	double mean_prior[k][d];
	double invwish_prior[d*d];
	double nu_prior;
	double wishlndet_prior;
	double lnC_prior;
	
	alpha_prior = pow(10.0, -2.0);
	beta_prior = 1.0;
	nu_prior = d;
	
	GetRNGstate();
	for(i=0; i<k; i++) {
		for(j=0; j<d; j++) {
			interm1 = 0.5 * (maxs[j] - mins[j]);
			mean_prior[i][j] = runif(mins[j] - interm1, maxs[j] + interm1);
		}
	}
	PutRNGstate();
	
	
	for(i=0; i<(d*d); i++) {
		invwish_prior[i] = 0.0;
		tempmat1[i] = 0.0;
	}
	

	//interm1 = d / REAL(coerceVector(var, REALSXP))[0]; var defaults to 128
	interm1 = d / 128.0;
	matview1 = gsl_matrix_view_array(invwish_prior, d, d);
	matview2 = gsl_matrix_view_array(tempmat1, d, d);
	for(i=0; i<d; i++) {
		interm2 = pow(maxs[i] - mins[i], 2.0) * interm1;
		gsl_matrix_set(&matview1.matrix, i, i, interm2);
		interm2 = 1.0/interm2;
		gsl_matrix_set(&matview2.matrix, i, i, interm2);
	}
	
	gsl_linalg_LU_decomp(&matview2.matrix, perm, &ind);
	wishlndet_prior = gsl_linalg_LU_lndet(&matview2.matrix);
	
	lnC_prior = 0.0;
	lnC_prior += gsl_sf_lngamma(k * alpha_prior);
	for(i=0; i<k; i++) {
		lnC_prior -= gsl_sf_lngamma(alpha_prior);
	}
		
	// declare estimated model
	double alpha_model[k];
	double beta_model[k];
	double mean_model[k][d];
	double wish_model[k][d*d];
	double wishlndet_model[k];
	double nu_model[k];
	double lnC_model;
	
	// init with prior values
	lnC_model = lnC_prior;
	for(i=0; i<k; i++) {
		for(j=0; j<(d*d); j++) {
			wish_model[i][j] = 0.0;
		}
	}
	matview1 = gsl_matrix_view_array(invwish_prior, d, d);
	for(i=0; i<k; i++) {
		alpha_model[i] = alpha_prior;
		beta_model[i] = beta_prior;
		nu_model[i] = nu_prior;
		wishlndet_model[i] = wishlndet_prior;
		matview2 = gsl_matrix_view_array(wish_model[i], d, d);
		for(j=0; j<d; j++) {
			mean_model[i][j] = mean_prior[i][j];
			gsl_matrix_set(&matview2.matrix, j, j, 1.0 / gsl_matrix_get(&matview1.matrix, j, j));
		}
	}
	
	// declare moments and statistics
	double lnPi[k];
	double lnLambda[k];
	
	double Nk[k];
	double meank[k][d];
	double Sk[k][d*d];
	double Ck[k][d*d];
	
	// calculation variables
	double lnresp[k];
	double resp[k];
	
	// start algo (current bound = GSL_NEGINF)
	double bound = GSL_NEGINF;
	double newbound;
	double pXlog = 0.0;
	double pZlog = 0.0;
	double pPiLog = 0.0;
	double pMuLog = 0.0;
	double qZlog = 0.0;
	double qPiLog = 0.0;
	double qMuLog = 0.0;
	
	// compute moments pi_k and lambda_k before iterations
	view1 = gsl_vector_view_array(alpha_model, k);
	//interm1 = sum(&view1.vector);
	interm1 = sum(alpha_model, k, 1);

	interm1 = gsl_sf_psi(interm1);
	for(i=0; i<k; i++) {
		lnPi[i] = gsl_sf_psi(alpha_model[i]) - interm1;
	}
	
	for(i=0; i<k; i++) {
		lnLambda[i] = d * gsl_sf_log(2.0) + wishlndet_model[i];
		for(j=0; j<d; j++) {
			lnLambda[i] += gsl_sf_psi((nu_model[i] - j) / 2.0);
		}
	}	
	

	unsigned converge = 0;
	unsigned itcount=1;
	
	Rprintf("init\n");
	
	while(!converge) {
			
		// first reinit statistics (as we will incrementally update them)
		for(i=0; i<k; i++) {
			Nk[i] = 0.0;
			for(j=0; j<d; j++) {
				meank[i][j] = 0.0;
				for(u=0; u<d; u++) {
					Sk[i][j*d + u] = 0.0;
					Ck[i][j*d + u] = 0.0;
				}
			}
		}
		
		// prepare pZlog and qZlog
		qZlog = 0.0;
		
		// for each l element :  for each k : compute the specific moment, the log resp.
		// incremental update to statistics		
		for(i=0; i<L; i++) {
			for(j=0; j<k; j++) {
				// first add up already know elements
				lnresp[j] = 2.0 * lnPi[j] + lnLambda[j] - d * gsl_sf_log(2.0 * M_PI);
				// trace element : matrix mult.
				matview1 = gsl_matrix_view_array(wish_model[j], d, d);
				matview2 = gsl_matrix_view_array(REAL(VECTOR_ELT(cov, i)), d, d);
				matview3 = gsl_matrix_view_array(tempmat1, d, d);
				gsl_blas_dsymm(CblasLeft, CblasUpper, nu_model[j], &matview1.matrix, &matview2.matrix, 0.0, &matview3.matrix);
				view1 = gsl_matrix_diagonal(&matview3.matrix);
				//lnresp[j] -= sum(&view1.vector);
				read = gsl_vector_ptr(&view1.vector, 0);
				lnresp[j] -= sum(read, d, view1.vector.stride);
				
				// quadratic term
				view1 = gsl_vector_view_array(tempvec1, d);
				view2 = gsl_vector_view_array(REAL(VECTOR_ELT(mean, i)), d);
				gsl_blas_dcopy(&view2.vector, &view1.vector);
				view2 = gsl_vector_view_array(mean_model[j], d);
				gsl_blas_daxpy(-1.0, &view2.vector, &view1.vector);
				matview1 = gsl_matrix_view_array(wish_model[j], d, d);
				view2 = gsl_vector_view_array(tempvec2, d);
				gsl_blas_dsymv(CblasUpper, nu_model[j], &matview1.matrix, &view1.vector, 0.0, &view2.vector);
				gsl_blas_ddot(&view1.vector, &view2.vector, &interm1);
				lnresp[j] -= interm1;
				lnresp[j] -= d / beta_model[j];
				
				lnresp[j] = lnresp[j] * (N * REAL(w)[i] / 2.0);
				
				// manage limits of gsl_sf_exp
				/* 2/03 : change norm strategy
				if(lnresp[j] < -700.0) {
					resp[j] = 0.0;
				} else {
					resp[j] = gsl_sf_exp(lnresp[j]);
				}*/
			}
			// when a line is finished, normalize responsibilities and update stats
			// 2/03 : adjust lnresp and compute exp in any case.
			interm1 = max_val(lnresp, k);
			
			for(j=0; j<k; j++) {
				lnresp[j] -= interm1;
				if(lnresp[j] > -700.0) {
					resp[j] = gsl_sf_exp(lnresp[j]);
				} else {
					resp[j] = 0.0;
				}
			}
			
			view1 = gsl_vector_view_array(resp, k);
			//interm1 = sum(&view1.vector);
			interm1 = sum(resp, k, 1);
			
			// if sum = 0, normalization is not possible : we resort to arg max lnresp.
			// else no problem.
			// 2/03 : useless now.
			/*if(interm1 == 0.0) {
				view2 = gsl_vector_view_array(lnresp, k);
				ind = gsl_vector_max_index(&view2.vector);
				for(j=0; j<k; j++) {
					if(j==ind) {
						resp[j] = 1.0;
					} else {
						resp[j] = 0.0;
					}
				}
			} else {
				*/
			for(j=0; j<k; j++) {
				resp[j] /= interm1;
			}
			
			
			
			// update statisticspXlog
			for(j=0; j<k; j++) {
				interm1 = N * REAL(w)[i] * resp[j];
				Nk[j] += interm1;
				view1 = gsl_vector_view_array(REAL(VECTOR_ELT(mean, i)), d);
				view2 = gsl_vector_view_array(meank[j], d);
				gsl_blas_daxpy(interm1, &view1.vector, &view2.vector);
				
				matview1 = gsl_matrix_view_array(Sk[j], d, d);
				gsl_blas_dsyr(CblasUpper, interm1, &view1.vector, &matview1.matrix);
				
				matview1 = gsl_matrix_view_array(Ck[j], d, d);
				matview2 = gsl_matrix_view_array(tempmat1, d, d);
				matview3 = gsl_matrix_view_array(REAL(VECTOR_ELT(cov, i)), d, d);
				
				gsl_matrix_memcpy(&matview2.matrix, &matview3.matrix);
				gsl_matrix_scale(&matview2.matrix, interm1);
				gsl_matrix_add(&matview1.matrix, &matview2.matrix);
			}
			
			
			// directly update qZlog as here resp are needed
			for(j=0; j<k; j++) {
				// manage case where resp[j] = 0 : continuity prolongation
				if(resp[j] != 0.0) {
					qZlog += resp[j] * gsl_sf_log(resp[j]);
				}
			}
			
		}
		
		// when finished :
		// adjust statistics for each k
		for(i=0; i<k; i++) {
			// scale meank by the final Nk
			// manage case where Nk is 0
			if(Nk[i] != 0.0) {
				view1 = gsl_vector_view_array(meank[i], d);
				gsl_blas_dscal(1.0 / Nk[i], &view1.vector);
				
				// update Sk : newSk = (1/Nk) * currentSk - meank.meank^T
				matview1 = gsl_matrix_view_array(Sk[i], d, d);
				gsl_matrix_scale(&matview1.matrix, 1.0 / Nk[i]);
				gsl_blas_dsyr(CblasUpper, -1.0, &view1.vector, &matview1.matrix);
				
				// update Ck
				matview1 = gsl_matrix_view_array(Ck[i], d, d);
				gsl_matrix_scale(&matview1.matrix, 1.0 / Nk[i]);
			} else {
				for(j=0; j<d; j++) {
					meank[i][j] = 0.0;
					for(u=0; u<d; u++) {
						Sk[i][j*d + u] = 0.0;
						Ck[i][j*d + u] = 0.0;
					}
				}
			}
		}
		
		// reinit model and bound before update
		for(i=0; i<k; i++) {
			alpha_model[i] = 0.0;
			beta_model[i] = 0.0;
			wishlndet_model[i] = 0.0;
			nu_model[i] = 0.0;
			for(j=0; j<d; j++) {
				mean_model[i][j] = 0.0;
				for(u=0; u<d; u++) {
					wish_model[i][j*d + u] = 0.0;
				}
			}
		}
		lnC_model = 0.0;
		pXlog = 0.0;
		pZlog = 0.0;
		pPiLog = 0.0;
		pMuLog = 0.0;
		qPiLog = 0.0;
		qMuLog = 0.0;
		

					
		// update model for each k,
		// then bound
		
		for(i=0; i<k; i++) {
			alpha_model[i] = alpha_prior + Nk[i];
			beta_model[i] = beta_prior + Nk[i];
			nu_model[i] = nu_prior + Nk[i];
		}
		// precompute element for lnPi
		view1 = gsl_vector_view_array(alpha_model, k);
		//interm1 = sum(&view1.vector);
		interm1 = sum(alpha_model, k, 1);
		interm1 = gsl_sf_psi(interm1);
		
		for(i=0; i<k; i++) {
			lnPi[i] = -interm1;
		}	
		
		for(i=0; i<k; i++) {
	
			// mean update
			view1 = gsl_vector_view_array(mean_model[i], d);
			view2 = gsl_vector_view_array(mean_prior[i], d);
			gsl_blas_daxpy(beta_prior, &view2.vector, &view1.vector);
			view2 = gsl_vector_view_array(meank[i], d);
			gsl_blas_daxpy(Nk[i], &view2.vector, &view1.vector);
			gsl_blas_dscal(1.0 / beta_model[i], &view1.vector);
		
			// inverse Wish update 1st term
			matview1 = gsl_matrix_view_array(wish_model[i], d, d);
			matview2 = gsl_matrix_view_array(invwish_prior, d, d);
			gsl_matrix_memcpy(&matview1.matrix, &matview2.matrix);
			
			// 2nd term
			matview2 = gsl_matrix_view_array(Sk[i], d, d);
			matview3 = gsl_matrix_view_array(tempmat1, d, d);
			gsl_matrix_memcpy(&matview3.matrix, &matview2.matrix);
			gsl_matrix_scale(&matview3.matrix, Nk[i]);
			gsl_matrix_add(&matview1.matrix, &matview3.matrix);
			
			// 3rd term
			view1 = gsl_vector_view_array(tempvec1, d);
			view2 = gsl_vector_view_array(meank[i], d);
			gsl_blas_dcopy(&view2.vector, &view1.vector);
			view2 = gsl_vector_view_array(mean_prior[i], d);
			gsl_blas_daxpy(-1.0, &view2.vector, &view1.vector);
			gsl_blas_dsyr(CblasUpper, (beta_prior * Nk[i]) / beta_model[i], &view1.vector, &matview1.matrix);
			
			// 4th term
			matview2 = gsl_matrix_view_array(tempmat1, d, d);
			matview3 = gsl_matrix_view_array(Ck[i], d, d);
			gsl_matrix_memcpy(&matview2.matrix, &matview3.matrix);
			gsl_matrix_scale(&matview2.matrix, Nk[i]);
			gsl_matrix_add(&matview1.matrix, &matview2.matrix);
			
			// upper complete before inverse, computing lndet and bound updates
			//PBR
			//matview2 = gsl_matrix_view_array(wish_model[i]);
			//upperComplete(wish_model[i], d);
			upperComplete(&matview1.matrix);
			matview2 = gsl_matrix_view_array(tempmat1, d, d);
			gsl_matrix_memcpy(&matview2.matrix, &matview1.matrix);
			gsl_linalg_LU_decomp(&matview2.matrix, perm, &ind);
			gsl_linalg_LU_invert(&matview2.matrix, perm, &matview1.matrix);
			gsl_matrix_memcpy(&matview2.matrix, &matview1.matrix);
			gsl_linalg_LU_decomp(&matview2.matrix, perm, &ind);
			wishlndet_model[i] = gsl_linalg_LU_lndet(&matview2.matrix);
			
			// compute moments pi_k and lambda_k before bound update
			lnPi[i] += gsl_sf_psi(alpha_model[i]);
		
			lnLambda[i] = d * gsl_sf_log(2.0) + wishlndet_model[i];
			for(j=0; j<d; j++) {
				lnLambda[i] += gsl_sf_psi((nu_model[i] - j) / 2.0);
			}			
			
			
			// bound updates
			// pXlog
			interm1 = lnLambda[i] - d / beta_model[i] - d * gsl_sf_log(2.0 * M_PI);
			matview2 = gsl_matrix_view_array(Sk[i], d, d);
			matview3 = gsl_matrix_view_array(tempmat1, d, d);
			gsl_matrix_memcpy(&matview3.matrix, &matview2.matrix);
			matview2 = gsl_matrix_view_array(Ck[i], d, d);
			gsl_matrix_add(&matview3.matrix, &matview2.matrix);
			matview2 = gsl_matrix_view_array(tempmat2, d, d);
			gsl_blas_dsymm(CblasLeft, CblasUpper, nu_model[i], &matview3.matrix, &matview1.matrix, 0.0, &matview2.matrix);
			view1 = gsl_matrix_diagonal(&matview2.matrix);
			read = gsl_vector_ptr(&view1.vector, 0);
			//interm1 -= sum(&view1.vector);
			interm1 -= sum(read, d, view1.vector.stride);
			
			view1 = gsl_vector_view_array(tempvec1, d);
			view2 = gsl_vector_view_array(meank[i], d);
			gsl_blas_dcopy(&view2.vector, &view1.vector);
			view2 = gsl_vector_view_array(mean_model[i], d);
			gsl_blas_daxpy(-1.0, &view2.vector, &view1.vector);
			view2 = gsl_vector_view_array(tempvec2, d);
			gsl_blas_dsymv(CblasUpper, nu_model[i], &matview1.matrix, &view1.vector, 0.0, &view2.vector);
			gsl_blas_ddot(&view1.vector, &view2.vector, &interm2);
			interm1 -= interm2;
			
			pXlog += 0.5 * Nk[i] * interm1;
			
			// pZlog
			pZlog += Nk[i] * lnPi[i];
			
			
			// pPiLog
			pPiLog += (alpha_prior - 1.0) * lnPi[i];
			
			// pMuLog
			interm1 = d * gsl_sf_log(beta_prior / (2.0 * M_PI)) + lnLambda[i] - (d * beta_prior) / beta_model[i];
			view1 = gsl_vector_view_array(tempvec1, d);
			view2 = gsl_vector_view_array(mean_model[i], d);
			gsl_blas_dcopy(&view2.vector, &view1.vector);
			view2 = gsl_vector_view_array(mean_prior[i], d);
			gsl_blas_daxpy(-1.0, &view2.vector, &view1.vector);
			view2 = gsl_vector_view_array(tempvec2, d);
			gsl_blas_dsymv(CblasUpper, beta_prior * nu_model[i], &matview1.matrix, &view1.vector, 0.0, &view2.vector);
			gsl_blas_ddot(&view1.vector, &view2.vector, &interm2);
			interm1 -= interm2;
			pMuLog += 0.5 * interm1;
			
			pMuLog += (nu_prior - d - 1) * lnLambda[i] / 2.0;
			
			matview2 = gsl_matrix_view_array(invwish_prior, d, d);
			matview3 = gsl_matrix_view_array(tempmat1, d, d);
			gsl_blas_dsymm(CblasLeft, CblasUpper, nu_model[i], &matview2.matrix, &matview1.matrix, 0.0, &matview3.matrix);
			view1 = gsl_matrix_diagonal(&matview3.matrix);
			read = gsl_vector_ptr(&view1.vector, 0);
			//interm1 = sum(&view1.vector);
			interm1 = sum(read, d, view1.vector.stride);
			pMuLog -= 0.5 * interm1;
			
			// qPiLog
			qPiLog += (alpha_model[i] - 1.0) * lnPi[i];
			
			// qMuLog
			qMuLog += 0.5 * lnLambda[i];
			qMuLog += (d * (gsl_sf_log(beta_model[i] / (2 * M_PI)) - 1.0)) / 2.0;
			qMuLog += -nu_model[i] * wishlndet_model[i] / 2.0 - nu_model[i] * d * gsl_sf_log(2.0) / 2.0 - d * (d-1) * gsl_sf_log(M_PI) / 4.0;
			for(j=0; j<d; j++) {
				qMuLog -= gsl_sf_lngamma((nu_model[i] - j) / 2.0);
			}
			qMuLog += ((nu_model[i] - d - 1.0) / 2.0) * lnLambda[i];
			qMuLog -= nu_model[i] * d / 2.0;
			
		}
		
		// final updates of the model (ln C alpha, lnCprior, pMuLog) outside of a K loop
		// and finalize pPiLog
		
		// calculate lnC_model
		view1 = gsl_vector_view_array(alpha_model, k);
		//interm1 = sum(&view1.vector);
		interm1 = sum(alpha_model, k, 1);
		lnC_model = gsl_sf_lngamma(interm1);
		for(i=0; i<k; i++) {
			lnC_model -= gsl_sf_lngamma(alpha_model[i]);
		}
		
		// bound final terms
		interm1 = -nu_prior * wishlndet_prior / 2.0 - nu_prior * d * gsl_sf_log(2.0) / 2.0 - d * (d-1) * gsl_sf_log(M_PI) / 4.0;
		for(j=0; j<d; j++) {
			interm1 -= gsl_sf_lngamma((nu_prior - j) / 2.0);
		}
		pMuLog += k * interm1;
		
		pPiLog += lnC_prior;
		
		qPiLog += lnC_model;
		
		// check convergence
		
		newbound = pXlog + pZlog + pPiLog + pMuLog - qZlog - qPiLog - qMuLog;
		
		if(maxit==R_NilValue) {
			if((newbound - bound) < REAL(coerceVector(thres, REALSXP))[0]) {
				converge = 1;
			}
		} else {
			if(itcount >= INTEGER(coerceVector(maxit, INTSXP))[0]) {
				converge = 1;
			}
		}

		bound = newbound;
		Rprintf("bound is %f\n", bound);
		
		itcount++;
	
	}
	
	
	// allocate and build returned data
	SEXP model_out, ret_out, alpha_out, beta_out, nu_out, mean_out, wish_out, curmean, curwish, dims1, dims2;
	PROTECT(alpha_out = allocVector(REALSXP, k));
	PROTECT(beta_out = allocVector(REALSXP, k));
	PROTECT(nu_out = allocVector(REALSXP, k));
	PROTECT(mean_out = allocVector(VECSXP, k));
	PROTECT(wish_out = allocVector(VECSXP, k));
	for(i=0; i<k; i++) {
		SET_VECTOR_ELT(mean_out, i, curmean = allocVector(REALSXP, d));
		SET_VECTOR_ELT(wish_out, i, curwish = allocMatrix(REALSXP, d, d));
		
		REAL(alpha_out)[i] = alpha_model[i];
		REAL(beta_out)[i] = beta_model[i];
		REAL(nu_out)[i] = nu_model[i];
		for(j=0; j<d; j++) {
			REAL(curmean)[j] = mean_model[i][j];
			for(u=0; u<d; u++) {
				REAL(curwish)[j*d + u] = wish_model[i][j*d + u];
			}
		}
	}
	
	SEXP nbits;
	PROTECT(nbits = allocVector(INTSXP, 1));
	INTEGER(nbits)[0] = itcount;
	
	PROTECT(model_out = allocVector(VECSXP, 6));
	SET_VECTOR_ELT(model_out, 0, alpha_out);
	SET_VECTOR_ELT(model_out, 1, beta_out);
	SET_VECTOR_ELT(model_out, 2, nu_out);
	SET_VECTOR_ELT(model_out, 3, mean_out);
	SET_VECTOR_ELT(model_out, 4, wish_out);
	SET_VECTOR_ELT(model_out, 5, nbits);
	
	PROTECT(ret_out = allocVector(VECSXP, 1));
	SET_VECTOR_ELT(ret_out, 0, model_out);
	
	
	PROTECT(dims1 = allocVector(VECSXP, 1));
	PROTECT(dims2 = allocVector(VECSXP, 6));
	SET_VECTOR_ELT(dims1, 0, mkChar("model"));
	SET_VECTOR_ELT(dims2, 0, mkChar("alpha"));
	SET_VECTOR_ELT(dims2, 1, mkChar("beta"));
	SET_VECTOR_ELT(dims2, 2, mkChar("nu"));
	SET_VECTOR_ELT(dims2, 3, mkChar("mean"));
	SET_VECTOR_ELT(dims2, 4, mkChar("wish"));
	SET_VECTOR_ELT(dims2, 5, mkChar("nbits"));
	
	setAttrib(model_out, R_NamesSymbol, dims2);
	setAttrib(ret_out, R_NamesSymbol, dims1);
	
		
	gsl_permutation_free(perm);
	UNPROTECT(10);
	return(ret_out);
		
}





