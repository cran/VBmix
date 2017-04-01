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


static double sumGSL(gsl_vector *vec) {
	// computes sum of a vector. Uses a pointer for quicker.
	int i;
	double tot = 0.0;
	double *cur = gsl_vector_ptr(vec, 0);
	for(i=0; i<vec->size; i++) {
		tot += cur[i*vec->stride];
	}
	
	return(tot);
}

static double sumDouble(double *vec, int size) {
	int i;
	double tot = 0.0;
	for(i=0; i<size; i++) {
		tot += vec[i];
	}
	return(tot);
}

static int max_index(double *vec, int size) {
	int i;
	double max=GSL_NEGINF;
	int maxind;
	for(i=0; i<size; i++) {
		if(vec[i] > max) {
			max = vec[i];
			maxind = i;
		}
	}
	return(maxind);
}

static int max_int_val(int *vec, int size) {
	int i;
	int themax = 0;
	for(i=0; i<size; i++) {
		themax = (vec[i] > themax) ? vec[i] : themax;
	}
	return(themax);
}

static double max_val(double *vec, int size) {
	int i;
	double themax = GSL_NEGINF;
	for(i=0; i<size; i++) {
		themax = (vec[i] > themax) ? vec[i] : themax;
	}
	return(themax);
}




// rho, var, ajouter maxit et thres
SEXP vbconstr(SEXP models, SEXP ncomp, SEXP thres, SEXP maxit, SEXP rho) {
	// for convenience get pointers to data structures
	// parameter a : L x P binary matrix denoting constraints between components in models
	SEXP w, mean, cov, a;
	w = getListElement(models, "w");
	mean = getListElement(models, "mean");
	cov = getListElement(models, "cov");
	a = getListElement(models, "a");
	// init prior model + precalc quantities (inverse wishart, constant terms in bound)
	// get sizes from data
	int L = length(w);
	int k = INTEGER(coerceVector(ncomp, INTSXP))[0];
	int d = length(VECTOR_ELT(mean, 0));
	int N = 200000;
	
	// for this to work, constraints have to be 0 0 ... n-1 n-1
	int nbim = max_int_val(INTEGER(coerceVector(a, INTSXP)), length(a)) + 1;
	
	// declare stack variables for all usages
	int i, j, u, v, ind;
	double interm1;
	double interm2;
	int *intptr1;
	int *intptr2;
	double *dbptr1;
	gsl_vector_view view1;

	gsl_matrix *tempmat1 = gsl_matrix_alloc(d,d);
	gsl_matrix *tempmat2 = gsl_matrix_alloc(d,d);
	gsl_vector *tempvec1 = gsl_vector_alloc(d);
	gsl_vector *tempvec2 = gsl_vector_alloc(d);
	gsl_permutation *perm = gsl_permutation_alloc(d);
	
	// copy input means and covs into gsl_vectors/matrices (=> shall be used in algebra only)
	// double pointers to the weight vectors.
	double *wPtr = REAL(w);
	gsl_vector **inMeans = Calloc(L, gsl_vector *);
	gsl_matrix **inCovs = Calloc(L, gsl_matrix *);
	for(i=0; i<L; i++) {
		inMeans[i] = gsl_vector_alloc(d);
		dbptr1 = gsl_vector_ptr(inMeans[i], 0);
		memcpy(dbptr1, REAL(VECTOR_ELT(mean, i)), d * sizeof(double));
		inCovs[i] = gsl_matrix_alloc(d,d);
		// R reads mats by columns and gsl by lines : no issue since 
		dbptr1 = gsl_matrix_ptr(inCovs[i], 0, 0);
		memcpy(dbptr1, REAL(VECTOR_ELT(cov, i)), d*d*sizeof(double));
	}
	
	
	// need to store mpk : create a structure
	double mpk[nbim][k];
	int mpkbis[nbim][k];

	// variables to implement constraints : 
	// random order
	// a converted to indices of maxline
	// variable that keeps track of max resp
	SEXP samp, labs;
	PROTECT(labs = allocVector(INTSXP, L));
	for(i=0; i<L; i++) {
		INTEGER(labs)[i] = i;
	}
	// We get directly the constraint vector in non matrix form
	int *constraints;
	int *order;
	//int maxresp[L];
	
	// sample to get order
	SEXP expr, it;

	PROTECT(it = expr = allocVector(LANGSXP, 3));
	
	SETCAR(it, install("sample"));
	it = CDR(it);
	SETCAR(it, labs);
	it = CDR(it);
	SETCAR(it, allocVector(INTSXP, 1));
	INTEGER(CAR(it))[0] = L;
	PROTECT(samp = eval(expr, rho));
	order = INTEGER(samp);
	
	// order now contains (for simple) 4 2 3 1 ...
	
	// constraints already in form 0 0 ... 0 1 1 ... etc
	constraints = INTEGER(coerceVector(a, INTSXP));
	
	
	// first get maxs and mins from the models
	double mins[d];
	double maxs[d];
	
	for(i=0; i<d; i++) {
		mins[i] = GSL_POSINF;
		maxs[i] = GSL_NEGINF;
	}
	
	// ici order pas important : stats de synthese
	for(i=0; i<L; i++) {
		for(j=0; j<d; j++) {
			if(gsl_vector_get(inMeans[i], j) > maxs[j]) {
				maxs[j] = gsl_vector_get(inMeans[i], j);
			}
			if(gsl_vector_get(inMeans[i], j) < mins[j]) {
				mins[j] = gsl_vector_get(inMeans[i], j);
			}
		}
	}
	
	// init estim model with prior model (+ ln det wishart + inv wish_0 + ln C alpha)
	// only prior means with k values : all other are identical for all k.
	double alpha_prior;
	double beta_prior;
	
	// mean and cov shall be used in algebra => GSL
	gsl_vector **mean_prior = Calloc(k, gsl_vector *);
	gsl_matrix *invwish_prior = gsl_matrix_alloc(d, d);
	double nu_prior;
	double wishlndet_prior;
	double lnC_prior;
	
	alpha_prior = pow(10.0, -2.0);
	beta_prior = 1.0;
	nu_prior = d;
	
	GetRNGstate();
	for(i=0; i<k; i++) {
		mean_prior[i] = gsl_vector_alloc(d);
		for(j=0; j<d; j++) {
			interm1 = 0.5 * (maxs[j] - mins[j]);
			gsl_vector_set(mean_prior[i], j, runif(mins[j] - interm1, maxs[j] + interm1));
		}
	}
	PutRNGstate();
	
	
	gsl_matrix_set_zero(invwish_prior);
	gsl_matrix_set_zero(tempmat1);
	
	interm1 = d / 128.0;
	for(i=0; i<d; i++) {
		interm2 = pow(maxs[i] - mins[i], 2.0) * interm1;
		gsl_matrix_set(invwish_prior, i, i, interm2);
		interm2 = 1.0/interm2;
		gsl_matrix_set(tempmat1, i, i, interm2);
	}
	
	gsl_matrix_memcpy(tempmat2, tempmat1);
	gsl_linalg_LU_decomp(tempmat1, perm, &ind);
	wishlndet_prior = gsl_linalg_LU_lndet(tempmat1);
	
	lnC_prior = 0.0;
	lnC_prior += gsl_sf_lngamma(k * alpha_prior);
	for(i=0; i<k; i++) {
		lnC_prior -= gsl_sf_lngamma(alpha_prior);
	}
		
	// declare estimated model
	double alpha_model[k];
	double beta_model[k];
	
	// mean and wish shall be used in algebra => GSL
	gsl_vector **mean_model = Calloc(k, gsl_vector *);
	gsl_matrix **wish_model = Calloc(k, gsl_matrix *);
	double wishlndet_model[k];
	double nu_model[k];
	double lnC_model;
	
	// init with prior values
	lnC_model = lnC_prior;

	for(i=0; i<k; i++) {
		alpha_model[i] = alpha_prior;
		beta_model[i] = beta_prior;
		nu_model[i] = nu_prior;
		wishlndet_model[i] = wishlndet_prior;
		// change for mean and wish
		mean_model[i] = gsl_vector_alloc(d);
		wish_model[i] = gsl_matrix_alloc(d,d);
		gsl_blas_dcopy(mean_prior[i], mean_model[i]);
		gsl_matrix_memcpy(wish_model[i], tempmat2);
	}
	
	// declare moments and statistics
	double lnPi[k];
	double lnLambda[k];
	
	double Nk[k];
	gsl_vector **meank = Calloc(k, gsl_vector *);
	gsl_matrix **Sk = Calloc(k, gsl_matrix *);
	gsl_matrix **Ck = Calloc(k, gsl_matrix *);
	for(i=0; i<k; i++) {
		meank[i] = gsl_vector_alloc(d);
		Sk[i] = gsl_matrix_alloc(d,d);
		Ck[i] = gsl_matrix_alloc(d,d);
	}
	
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
	// new : term for the constraints distribution (to integrate in the lower bound)
	double pazlog = 0.0;
	
	// compute moments pi_k and lambda_k before iterations
	interm1 = sumDouble(alpha_model, k);
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
	
	Rprintf("init\n");
	
	// init count of number of iterations
	unsigned itcount = 1;
	
	while(!converge) {
			
		// first reinit statistics (as we will incrementally update them)
		for(i=0; i<k; i++) {
			Nk[i] = 0.0;
			gsl_vector_set_zero(meank[i]);
			gsl_matrix_set_zero(Sk[i]);
			gsl_matrix_set_zero(Ck[i]);
			for(j=0; j<nbim; j++) {
				mpk[j][i] = 0.0;
				mpkbis[j][i] = 0;
			}
		}
		
		// prepare qZlog
		qZlog = 0.0;
		
		// for each l element :  for each k : compute the specific moment, the log resp.
		// incremental update to statistics		
		for(i=0; i<L; i++) {
			for(j=0; j<k; j++) {

				// first add up already know elements
				lnresp[j] = 2.0 * lnPi[j] + lnLambda[j] - d * gsl_sf_log(2.0 * M_PI);
				// trace element : matrix mult.
				gsl_blas_dsymm(CblasLeft, CblasUpper, nu_model[j], wish_model[j], inCovs[order[i]], 0.0, tempmat1);
				view1 = gsl_matrix_diagonal(tempmat1);
				lnresp[j] -= sumGSL(&view1.vector);
				
				// quadratic term
				gsl_blas_dcopy(inMeans[order[i]], tempvec1);
				gsl_blas_daxpy(-1.0, mean_model[j], tempvec1);
				gsl_blas_dsymv(CblasUpper, nu_model[j], wish_model[j], tempvec1, 0.0, tempvec2);
				gsl_blas_ddot(tempvec1, tempvec2, &interm1);
				lnresp[j] -= interm1;
				lnresp[j] -= d / beta_model[j];
				
				lnresp[j] = lnresp[j] * (N * wPtr[order[i]] / 2.0);
				
				// additonnal section : we manage constraints
				/* change strategy : prevent explosion by using mpk (or mpkbis, discrete version) elegantly
				interm2 = 1.0;				
				for(u=0; u<i; u++) {
					if(maxresp[u] == j) {
						if(constraints[order[u]] == constraints[order[i]]) {
							interm2 += i-u+1;
						}
					}
				}*/
				ind = mpkbis[constraints[order[i]]][j] + 1;
				interm2 = ind * (ind + 1.0) / 2.0;
				mpkbis[constraints[order[i]]][j]++;
				
				
				// add constraints to current estimate
				// ress avec pow = 3
				// ress2 avec pow=2
				// ress3 avec pow=1
				lnresp[j] -= pow(interm2, 2.0);
				
				
				/* manage limits of gsl_sf_exp
				if(lnresp[j] < -700.0) {
					resp[j] = 0.0;
				} else if(lnresp[j] > 700.0) {
					resp[j] = GSL_POSINF;
				} else {
					resp[j] = gsl_sf_exp(lnresp[j]);
				}*/
			}
			
			// 2/03 : renormalize lnresp
			interm1 = max_val(lnresp, k);
			for(j=0; j<k; j++) {
				lnresp[j] -= interm1;
				if(lnresp[j] < -700.0) {
					resp[j] = 0.0;
				} else {
					resp[j] = gsl_sf_exp(lnresp[j]);
				}
			}
			
			
			
			// when a line is finished, normalize responsibilities and update stats
			interm1 = sumDouble(resp, k);
			interm2 = GSL_NEGINF;
			/* 2/03 : useless now
			// if sum = 0, normalization is not possible : we resort to arg max lnresp.
			// else no problem.
			
			// addtional patch : manage maxresp

			if((interm1 == 0.0) || (interm1 == GSL_POSINF)) {
				ind = max_index(lnresp, k);
				for(j=0; j<k; j++) {
					if(j==ind) {
						resp[j] = 1.0;
						maxresp[i] = j;
					} else {
						resp[j] = 0.0;
					}
				}
			} else {
				for(j=0; j<k; j++) {
					resp[j] /= interm1;
					if(lnresp[j] > interm2) {
						interm2 = lnresp[j];
						maxresp[i] = j;
					}
				}
			}
			*/
			
			for(j=0; j<k; j++) {
				resp[j] /= interm1;
				if(lnresp[j] > interm2) {
					interm2 = lnresp[j];
					// not used anymore
					//maxresp[i] = j;
				}
			}			
			
			// update mpk
			for(j=0; j<k; j++) {
				for(u=0; u<nbim; u++) {
					if(constraints[order[i]] == u) {
						mpk[u][j] += resp[j];
					}
				}
			}
			
			// update statistics
			for(j=0; j<k; j++) {
				interm1 = N * wPtr[order[i]] * resp[j];
				Nk[j] += interm1;
				gsl_blas_daxpy(interm1, inMeans[order[i]], meank[j]);
				
				gsl_blas_dsyr(CblasUpper, interm1, inMeans[order[i]], Sk[j]);
				
				gsl_matrix_memcpy(tempmat1, inCovs[order[i]]);
				gsl_matrix_scale(tempmat1, interm1);
				gsl_matrix_add(Ck[j], tempmat1);
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
				gsl_blas_dscal(1.0 / Nk[i], meank[i]);
				
				// update Sk : newSk = (1/Nk) * currentSk - meank.meank^T
				gsl_matrix_scale(Sk[i], 1.0 / Nk[i]);
				gsl_blas_dsyr(CblasUpper, -1.0, meank[i], Sk[i]);
				
				// update Ck
				gsl_matrix_scale(Ck[i], 1.0 / Nk[i]);
			}
			/* useless 
			else {
				for(j=0; j<d; j++) {
					meank[i][j] = 0.0;
					for(u=0; u<d; u++) {
						Sk[i][j*d + u] = 0.0;
						Ck[i][j*d + u] = 0.0;
					}
				}
			}*/
		}
		
		// reinit model and bound before update
		for(i=0; i<k; i++) {
			alpha_model[i] = 0.0;
			beta_model[i] = 0.0;
			wishlndet_model[i] = 0.0;
			nu_model[i] = 0.0;
			gsl_vector_set_zero(mean_model[i]);
			gsl_matrix_set_zero(wish_model[i]);
		}
		lnC_model = 0.0;
		pXlog = 0.0;
		pZlog = 0.0;
		pPiLog = 0.0;
		pMuLog = 0.0;
		qPiLog = 0.0;
		qMuLog = 0.0;
		pazlog = - k * nbim;
		

					
		// update model for each k,
		// then bound
		
		for(i=0; i<k; i++) {
			alpha_model[i] = alpha_prior + Nk[i];
			beta_model[i] = beta_prior + Nk[i];
			nu_model[i] = nu_prior + Nk[i];
		}
		// precompute element for lnPi
		interm1 = sumDouble(alpha_model, k);
		interm1 = gsl_sf_psi(interm1);
		
		for(i=0; i<k; i++) {
			lnPi[i] = -interm1;
		}	
		
		for(i=0; i<k; i++) {
	
			// mean update
			gsl_blas_daxpy(beta_prior, mean_prior[i], mean_model[i]);
			gsl_blas_daxpy(Nk[i], meank[i], mean_model[i]);
			gsl_blas_dscal(1.0 / beta_model[i], mean_model[i]);
		
			// inverse Wish update 1st term
			gsl_matrix_memcpy(wish_model[i], invwish_prior);
			
			// 2nd term
			gsl_matrix_memcpy(tempmat1, Sk[i]);
			gsl_matrix_scale(tempmat1, Nk[i]);
			gsl_matrix_add(wish_model[i], tempmat1);
			
			// 3rd term
			gsl_blas_dcopy(meank[i], tempvec1);
			gsl_blas_daxpy(-1.0, mean_prior[i], tempvec1);
			gsl_blas_dsyr(CblasUpper, (beta_prior * Nk[i]) / beta_model[i], tempvec1, wish_model[i]);
			
			// 4th term
			gsl_matrix_memcpy(tempmat1, Ck[i]);
			gsl_matrix_scale(tempmat1, Nk[i]);
			gsl_matrix_add(wish_model[i], tempmat1);
			
			// upper complete before inverse, computing lndet and bound updates
			upperComplete(wish_model[i]);
			gsl_matrix_memcpy(tempmat1, wish_model[i]);
			// det(A^-1) = 1/det(A)
			// => ln det(A^-1) = - ln det(A)
			
			gsl_linalg_LU_decomp(tempmat1, perm, &ind);
			gsl_linalg_LU_invert(tempmat1, perm, wish_model[i]);
			wishlndet_model[i] = - gsl_linalg_LU_lndet(tempmat1);
			
			
			// compute moments pi_k and lambda_k before bound update
			lnPi[i] += gsl_sf_psi(alpha_model[i]);
		
			lnLambda[i] = d * gsl_sf_log(2.0) + wishlndet_model[i];
			for(j=0; j<d; j++) {
				lnLambda[i] += gsl_sf_psi((nu_model[i] - j) / 2.0);
			}			
			
			
			// bound updates
			// pXlog
			interm1 = lnLambda[i] - d / beta_model[i] - d * gsl_sf_log(2.0 * M_PI);
			gsl_matrix_memcpy(tempmat1, Sk[i]);
			gsl_matrix_add(tempmat1, Ck[i]);
			
			
			gsl_blas_dsymm(CblasLeft, CblasUpper, nu_model[i], tempmat1, wish_model[i], 0.0, tempmat2);
			view1 = gsl_matrix_diagonal(tempmat2);
			interm1 -= sumGSL(&view1.vector);
			
			gsl_blas_dcopy(meank[i], tempvec1);
			gsl_blas_daxpy(-1.0, mean_model[i], tempvec1);
			gsl_blas_dsymv(CblasUpper, nu_model[i], wish_model[i], tempvec1, 0.0, tempvec2);
			gsl_blas_ddot(tempvec1, tempvec2, &interm2);
			interm1 -= interm2;
			
			
			pXlog += 0.5 * Nk[i] * interm1;
			
			// pZlog
			pZlog += Nk[i] * lnPi[i];
			
			// pazlog
			for(j=0; j<nbim; j++) {
				v = (int)round(mpk[j][i]);
				for(u=0; u<=v; u++) {
					pazlog -= gsl_sf_log(1.0 + (double)u);
				}
			}
			
			// pPiLog
			pPiLog += (alpha_prior - 1.0) * lnPi[i];
			
			// pMuLog
			interm1 = d * gsl_sf_log(beta_prior / (2.0 * M_PI)) + lnLambda[i] - (d * beta_prior) / beta_model[i];
			gsl_blas_dcopy(mean_model[i], tempvec1);
			gsl_blas_daxpy(-1.0, mean_prior[i], tempvec1);
			gsl_blas_dsymv(CblasUpper, beta_prior * nu_model[i], wish_model[i], tempvec1, 0.0, tempvec2);
			gsl_blas_ddot(tempvec1, tempvec2, &interm2);
			interm1 -= interm2;
			pMuLog += 0.5 * interm1;
			
			pMuLog += (nu_prior - d - 1) * lnLambda[i] / 2.0;
			
			gsl_blas_dsymm(CblasLeft, CblasUpper, nu_model[i], invwish_prior, wish_model[i], 0.0, tempmat1);
			view1 = gsl_matrix_diagonal(tempmat1);
			interm1 = sumGSL(&view1.vector);
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
		interm1 = sumDouble(alpha_model, k);
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
		// modif bound
		newbound = pXlog + pZlog + pazlog + pPiLog + pMuLog - qZlog - qPiLog - qMuLog;
		if(maxit==R_NilValue) {
			if(newbound - bound < REAL(coerceVector(thres, REALSXP))[0]) {
				converge = 1;
			}
		} else {
			if(itcount >= INTEGER(coerceVector(maxit, INTSXP))[0]) {
				converge = 1;
			}
		}
		bound = newbound;
		Rprintf("bound is %f\n", bound);
		
		// update itcount
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
			dbptr1 = gsl_vector_ptr(mean_model[i], 0);
			memcpy(REAL(curmean), dbptr1, d*sizeof(double));
			dbptr1 = gsl_matrix_ptr(wish_model[i], 0, 0);
			memcpy(REAL(curwish), dbptr1, d*d*sizeof(double));
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
	
	
	// deallocate heap memory
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
	gsl_matrix_free(tempmat1);
	gsl_matrix_free(tempmat2);
	gsl_vector_free(tempvec1);
	gsl_vector_free(tempvec2);
	for(i=0; i<L; i++) {
		gsl_vector_free(inMeans[i]);
		gsl_matrix_free(inCovs[i]);
	}
	Free(inMeans);
	Free(inCovs);
	gsl_matrix_free(invwish_prior);
	for(i=0; i<k; i++) {
		gsl_vector_free(mean_prior[i]);
		gsl_vector_free(mean_model[i]);
		gsl_matrix_free(wish_model[i]);
		gsl_vector_free(meank[i]);
		gsl_matrix_free(Sk[i]);
		gsl_matrix_free(Ck[i]);
	}
	Free(mean_prior);
	Free(mean_model);
	Free(wish_model);
	Free(meank);
	Free(Sk);
	Free(Ck);
	
	UNPROTECT(13);
	
	
	return(ret_out);
		
}





