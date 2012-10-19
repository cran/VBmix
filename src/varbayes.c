// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include "utils.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_log.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_linalg.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


static void printvec(gsl_vector *vec, int cell) {
	int i;
	for(i=0; i<cell; i++) {
		Rprintf("%f ", gsl_vector_get(vec, i));
	}
	Rprintf("\n");
}

static void printmat(gsl_matrix *mat, int row, int col) {
	int i,j;
	for(i=0; i<row; i++) {
		for(j=0; j<col; j++) {
			Rprintf("%e ", gsl_matrix_get(mat, i, j));
		}
		Rprintf("\n");
	}
}





static void uppercomplete(gsl_matrix *mat) {
	int i,j;
	int n = mat->size1;
	for(i=0; i<(n-1); i++) {
		for(j=i+1; j<n; j++) {
			gsl_matrix_set(mat, j, i, gsl_matrix_get(mat, i, j));
		}
	}
}






	

static double sd(double *tab, double mean, int size, int stride) {
	int i;
	double tot=0.0;
	for(i=0; i<size; i++) {
		double interm = pow(tab[i*stride] - mean, 2);
		tot += interm;
	}
	tot = tot/(size - 1);
	tot = sqrt(tot);
	return(tot);
}


SEXP varbayes(SEXP data, SEXP ncomp, SEXP thres, SEXP maxit) {
	// 1 : coerce SEXP data into gsl_matrix
	
	PROTECT(data = coerceVector(data, REALSXP));
	PROTECT(ncomp = coerceVector(ncomp, INTSXP));
	
	double *read = REAL(data);
	
	double val;
	int i, j, m, ind;
	
	int n = INTEGER(getAttrib(data, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(data, R_DimSymbol))[1];
	int k = INTEGER(ncomp)[0];

	gsl_matrix *dat = gsl_matrix_alloc(n, d);
	
	for(i=0; i<n; i++) {
		for(j=0; j<d; j++) {
			gsl_matrix_set(dat, i, j, read[i + n*j]);
		}
	}

	double agitation[k];


	// 2 : extract basic stats in view of constructing prior
	double means[d];
	double sds[d];
	double mins[d];
	double maxs[d];
	
	double *tab;
	

	SEXP alpha;
	PROTECT(alpha=allocVector(REALSXP, k));
	for(i=0; i<k; i++) {
		REAL(alpha)[i] = pow(10, -2);
	}
	
	// threshold value
	double threshold = REAL(coerceVector(thres, REALSXP))[0];
	
	for(i=0; i<d; i++) {
		gsl_vector_view view1 = gsl_matrix_column(dat, i);
		double *val = gsl_vector_ptr(&view1.vector, 0);
		
		means[i] = sum(val, n, view1.vector.stride) / n;
		mins[i] = gsl_vector_min(&view1.vector);
		maxs[i] = gsl_vector_max(&view1.vector);
				
	}
	
	// 3 : build priors
	gsl_vector *alpha_prior = gsl_vector_alloc(k);
	gsl_vector *beta_prior = gsl_vector_alloc(k);
	gsl_vector *nu_prior = gsl_vector_alloc(k);
	gsl_vector **mean_prior = malloc(k * sizeof(gsl_vector *));
	gsl_matrix **wish_prior = malloc(k * sizeof(gsl_matrix *));
	for(i=0; i<k; i++) {
		mean_prior[i] = gsl_vector_alloc(d);
		wish_prior[i] = gsl_matrix_calloc(d,d);
	}
	
	for(i=0; i<k; i++) {
		double interm;
		gsl_vector_set(alpha_prior, i, REAL(alpha)[i]);
		gsl_vector_set(beta_prior, i, 1.0);
		gsl_vector_set(nu_prior, i, d);
		GetRNGstate();
		for(j=0; j<d; j++) {
			double samp = runif(mins[j], maxs[j]);
			gsl_vector_set(mean_prior[i], j, samp);
			// 12/08 : plus d'init selon variance
			interm = pow(maxs[j] - mins[j], 2);
			double tmp = interm * d;
			tmp = 64.0 / tmp;
			gsl_matrix_set(wish_prior[i], j, j, tmp);
		}
		PutRNGstate();
	}
	
	
	// intermediary : full allocation step => same objects used during all process
	gsl_matrix *c = gsl_matrix_calloc(n, d);
	gsl_matrix *betamatrix = gsl_matrix_alloc(n, k);
	gsl_vector *comput1 = gsl_vector_alloc(n);
	gsl_matrix *gammamat = gsl_matrix_alloc(k, d);
	gsl_matrix *decomp = gsl_matrix_alloc(d,d);
	gsl_matrix *pimatrix = gsl_matrix_alloc(n, k);
	gsl_matrix *kappamatrix = gsl_matrix_alloc(n, k);
	val = (d / 2) * gsl_sf_log(2 * M_PI);
	gsl_matrix *constmatrix = gsl_matrix_alloc(n, k);
	gsl_matrix_set_all(constmatrix, val);
	gsl_matrix *meanmatrix = gsl_matrix_alloc(d, n);
	gsl_vector *temp = gsl_vector_alloc(d);
	gsl_matrix *temp1 = gsl_matrix_calloc(d,d);
	gsl_matrix *temp2 = gsl_matrix_calloc(d,d);
	gsl_matrix *interm = gsl_matrix_calloc(d, n);
	gsl_matrix *prod = gsl_matrix_alloc(d,d);
	gsl_vector *meanvec = gsl_vector_alloc(d);
	gsl_matrix *calc = gsl_matrix_alloc(n, k);
	gsl_permutation *perm = gsl_permutation_alloc(d);
	gsl_matrix *inverse = gsl_matrix_alloc(d,d);
	gsl_vector *tempvec = gsl_vector_calloc(d);
	gsl_vector *templine = gsl_vector_calloc(k);
	
	
	
	
	// 4 : create and initiate moments object
	gsl_matrix *quad = gsl_matrix_calloc(n, k);
	gsl_vector *kappalog = gsl_vector_alloc(k);
	gsl_vector *pilog = gsl_vector_alloc(k);
	
	// 4.1 : compute quad	
	// new quad block : only 1 tempmat and meandat objects

	for(i=0; i<n; i++) {
		gsl_vector_view view1 = gsl_matrix_row(dat, i);
		for(j=0; j<k; j++) {
			gsl_vector_memcpy(meanvec, &view1.vector);
			gsl_vector_sub(meanvec, mean_prior[j]);
			gsl_blas_dsymv(CblasUpper, gsl_vector_get(nu_prior, j), wish_prior[j], meanvec, 0.0, temp);
			gsl_blas_ddot(temp, meanvec, &val);
			val += d / gsl_vector_get(beta_prior, j);
			gsl_matrix_set(quad, i, j, val);
		}
	}

	
	// 4.2 : compute kappalog
	for(i=0; i<d; i++) {
		gsl_vector_view view1 = gsl_matrix_column(gammamat, i);
		for(j=0; j<k; j++) {
			val = (gsl_vector_get(nu_prior, j) - i) / 2.0;
			gsl_vector_set(&view1.vector, j, digamma(val));
		}
	}
	
	
	
	

	
	for(i=0; i<k; i++) {
		gsl_vector_view view1 = gsl_matrix_row(gammamat, i);
		tab = gsl_vector_ptr(&view1.vector, 0);
		val = sum(tab, d, view1.vector.stride);
		val += d * gsl_sf_log(2.0);
		
		gsl_matrix_memcpy(decomp, wish_prior[i]);
		gsl_linalg_LU_decomp(decomp, perm, &ind);
		
		val += gsl_linalg_LU_lndet(decomp);
		gsl_vector_set(kappalog, i, val);
	}
	
	
	
	// 4.3 : compute pilog
	tab = gsl_vector_ptr(alpha_prior, 0);
	double pisum = sum(tab, k, alpha_prior->stride);
	
	pisum = digamma(pisum);
	for(i=0; i<k; i++) {
		val = digamma(gsl_vector_get(alpha_prior, i)) - pisum;
		gsl_vector_set(pilog, i, val);
	}
	
	
	// 5 computing statistics
	
	// 5.1 : building pimatrix and kappamatrix
	
	for(i=0; i<n; i++) {
		gsl_vector_view view1 = gsl_matrix_row(pimatrix, i);
		gsl_vector_memcpy(&view1.vector, pilog);
		view1 = gsl_matrix_row(kappamatrix, i);
		gsl_vector_memcpy(&view1.vector, kappalog);
	}
	
	
	// 5.2 building resp
	
	gsl_matrix *resp = gsl_matrix_calloc(n, k);
	gsl_matrix *oldresp = gsl_matrix_calloc(n, k);
	
	gsl_matrix_add(resp, pimatrix);

	gsl_matrix_scale(kappamatrix, 0.5);
	gsl_matrix_scale(quad, 0.5);
	
	gsl_matrix_add(resp, kappamatrix);
	gsl_matrix_sub(resp, constmatrix);
	gsl_matrix_sub(resp, quad);
	
	
	for(i=0; i<n; i++) {
		gsl_vector_view view1 = gsl_matrix_row(resp, i);
		tab = gsl_vector_ptr(&view1.vector, 0);
		
		// instead, detect highest value in lnresp, and normalize by these value before computing exp
		// => no special case any more
		val = gsl_vector_max(&view1.vector);
		for(j=0; j<k; j++) {
			gsl_vector_set(&view1.vector, j, exp(gsl_vector_get(&view1.vector, j) - val));
		}

		val = sum(tab, k, view1.vector.stride);
		for(j=0; j<k; j++) {
			tab[j * view1.vector.stride] = tab[j * view1.vector.stride] / val;
		}
		
		
		
	}
	

	
	// 5.3 building statistics
	
	gsl_vector *nk = gsl_vector_alloc(k);
	gsl_vector **meank = malloc(k * sizeof(gsl_vector *));
	gsl_matrix **sk = malloc(k * sizeof(gsl_matrix *));
	
	for(i=0; i<k; i++) {
		gsl_vector_view view1 = gsl_matrix_column(resp, i);
		tab = gsl_vector_ptr(&view1.vector, 0);
		val = sum(tab, n, view1.vector.stride);
		gsl_vector_set(nk, i, val);
		
		meank[i] = gsl_vector_calloc(d);
		// 12/08 : manage singularity
		if(val != 0.0) {
			gsl_blas_dgemv(CblasTrans, 1.0/val, dat, &view1.vector, 0.0, meank[i]);
		}
		
		
		sk[i] = gsl_matrix_calloc(d,d);
		
		if(val != 0.0) {
			// new (or old) strategy
			for(j=0; j<n; j++) {
				gsl_vector_view view2 = gsl_matrix_row(dat, j);
				gsl_vector_memcpy(meanvec, &view2.vector);
				gsl_vector_sub(meanvec, meank[i]);
				gsl_blas_dsyr(CblasUpper, gsl_vector_get(&view1.vector, j), meanvec, sk[i]);
			}
			
			gsl_matrix_scale(sk[i], 1.0/gsl_vector_get(nk, i));
			
			uppercomplete(sk[i]);
		} 
	
	}
		

	// 6 compute bound value
	
	// 6.1 pxlog
	
	double pxlog = 0.0;
	for(i=0; i<k; i++) {
		double tmp = 0.0;
		tmp += gsl_vector_get(kappalog, i);
		tmp -= d / gsl_vector_get(beta_prior, i);
		gsl_blas_dsymm(CblasRight, CblasUpper, gsl_vector_get(nu_prior, i), wish_prior[i], sk[i],  0.0, prod);
		gsl_vector_view view1 = gsl_matrix_diagonal(prod);
		tab = gsl_vector_ptr(&view1.vector, 0);
		val = sum(tab, d, view1.vector.stride);
		tmp -= val;
		gsl_vector_memcpy(meanvec, meank[i]);
		gsl_vector_sub(meanvec, mean_prior[i]);
		gsl_blas_dgemv(CblasNoTrans, gsl_vector_get(nu_prior, i), wish_prior[i], meanvec, 0.0, temp);
		gsl_vector_mul(temp, meanvec);
		tab = gsl_vector_ptr(temp, 0);
		val = sum(tab, d, temp->stride);
		tmp -= val;
		val = d * gsl_sf_log(2.0 * M_PI);
		tmp -= val;
		tmp = tmp * 0.5 * gsl_vector_get(nk, i);
		pxlog += tmp;
	}
		
	// 6.2 pzlog
	gsl_matrix_memcpy(calc, resp);
	gsl_matrix_mul_elements(calc, pimatrix);
	tab = gsl_matrix_ptr(calc, 0, 0);
	double pzlog = sum(tab, n*k, 1);
	
	
	// 6.3 ppilog
	double ppilog = 0.0;
	tab = gsl_vector_ptr(alpha_prior, 0);
	val = sum(tab, k, alpha_prior->stride);
	ppilog = gsl_sf_lngamma(val);
	for(i=0; i<k; i++) {
		ppilog -= gsl_sf_lngamma(gsl_vector_get(alpha_prior, i));
	}
	tab = gsl_vector_ptr(pilog, 0);
	val = sum(tab, k, pilog->stride);
	ppilog += (gsl_vector_get(alpha_prior, 0) - 1.0) * val;
	
	
	// 6.4 pmulog
	double pmulog = 0.0;
	for(i=0; i<k; i++) {
		double tmp = 0.0;
		tmp += d * gsl_sf_log(gsl_vector_get(beta_prior, 0) / (2.0 * M_PI));
		tmp += gsl_vector_get(kappalog, i);
		tmp -= d * gsl_vector_get(beta_prior, 0) / gsl_vector_get(beta_prior, i);
		// simplification : produit quadratique nul
		pmulog += 0.5 * tmp;
	}
	
	
	
	gsl_matrix_memcpy(decomp, wish_prior[0]);
	gsl_linalg_LU_decomp(decomp, perm, &ind);
		
	val = -(gsl_vector_get(nu_prior, 0) / 2.0) * gsl_linalg_LU_lndet(decomp);
	val -= gsl_vector_get(nu_prior, 0) * d * gsl_sf_log(2) / 2.0;
	val -= d * (d-1) * gsl_sf_log(M_PI) / 4.0;
	for(i=0; i<d; i++) {
		val -= gsl_sf_lngamma((gsl_vector_get(nu_prior, 0) - i) / 2.0);
	}
	
	pmulog += k * val;
	
	
	
	
	tab = gsl_vector_ptr(kappalog, 0);
	val = sum(tab, k, kappalog->stride);
	
	pmulog += (gsl_vector_get(nu_prior, 0) - d - 1.0) * val / 2.0;
	
	
	
	for(i=0; i<k; i++) {
		
		gsl_matrix_memcpy(decomp, wish_prior[i]);
		gsl_linalg_LU_decomp(decomp, perm, &ind);
		gsl_linalg_LU_invert(decomp, perm, inverse);
		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverse, wish_prior[i], 0.0, decomp);
		gsl_vector_view view1 = gsl_matrix_diagonal(decomp);
		tab = gsl_vector_ptr(&view1.vector, 0);
		
		pmulog -= 0.5 * gsl_vector_get(nu_prior, 0) * sum(tab, d, view1.vector.stride);
		
	}
	
	// 6.5 qzlog
	double qzlog = 0.0;
	tab = gsl_matrix_ptr(resp, 0, 0);
	for(i=0; i<=(n-1)*(k-1); i++) {
		if(tab[i] != 0) {
			qzlog += tab[i] * gsl_sf_log(tab[i]);
		}
	}
	
	// 6.6 qpilog
	double qpilog = 0.0;
	for(i=0; i<k; i++) {
		qpilog += (gsl_vector_get(alpha_prior, 0) - 1.0) * gsl_vector_get(pilog, i);
	}
	
	
	tab = gsl_vector_ptr(alpha_prior, 0);
	val = sum(tab, k, alpha_prior->stride);
	qpilog += gsl_sf_lngamma(val);
	for(i=0; i<k; i++) {
		qpilog -= gsl_sf_lngamma(gsl_vector_get(alpha_prior, i));
	}
	
	// 6.7 qmulog
	double qmulog = 0.0;

	
	for(i=0; i<k; i++) {
		double L = 0.0;
		double blog;
		for(j=0; j<d; j++) {
			L += gsl_sf_psi((gsl_vector_get(nu_prior, i) - j) / 2.0);
		}
		gsl_matrix_memcpy(decomp, wish_prior[i]);
		gsl_linalg_LU_decomp(decomp, perm, &ind);
		L += gsl_linalg_LU_lndet(decomp) + d * gsl_sf_log(2.0);
		
		blog = -(gsl_vector_get(nu_prior, 0) / 2.0) * gsl_linalg_LU_lndet(decomp);
		blog -= gsl_vector_get(nu_prior, 0) * d * gsl_sf_log(2) / 2.0;
		blog -= d * (d-1) * gsl_sf_log(M_PI) / 4.0;
		for(j=0; j<d; j++) {
			blog -= gsl_sf_lngamma((gsl_vector_get(nu_prior, 0) - j) / 2.0);
		}
		
		val = -((gsl_vector_get(nu_prior, i) - d - 1.0) / 2.0) * L + (gsl_vector_get(nu_prior, i) * d) / 2 - blog;
		
		qmulog += 0.5 * gsl_vector_get(kappalog, i) + d * gsl_sf_log(gsl_vector_get(beta_prior, i) / (2.0 * M_PI)) / 2.0 - d / 2.0 - val;
		
	}

	double bound = pxlog + pzlog + ppilog + pmulog - qzlog - qpilog - qmulog;
	double newbound;
	
	int count = 1;
	Rprintf("it %d bound is %f\n", count, bound);

	
	// create data for model to estimate
	gsl_vector *alpha_model = gsl_vector_alloc(k);
	gsl_vector *beta_model = gsl_vector_alloc(k);
	gsl_vector *nu_model = gsl_vector_alloc(k);
	gsl_vector **mean_model = malloc(k * sizeof(gsl_vector *));
	gsl_matrix **wish_model = malloc(k * sizeof(gsl_matrix *));
	for(i=0; i<k; i++) {
		mean_model[i] = gsl_vector_alloc(d);
		wish_model[i] = gsl_matrix_alloc(d,d);
	}

	int converge = 0;
	
	
	// build structures to store nk and agitation
	double **nks;
	double **agitations;

	
	// loop while converge is false
	while(!converge) {
		// 7 compute updated model
		for(i=0; i<k; i++) {
			double nki = gsl_vector_get(nk, i);
			double betai = gsl_vector_get(beta_prior, i);
			
			// 7.1 alpha, beta and nu
			gsl_vector_set(alpha_model, i, gsl_vector_get(alpha_prior, i) + nki);
			gsl_vector_set(beta_model, i, betai + nki);
			gsl_vector_set(nu_model, i, gsl_vector_get(nu_prior, i) + nki);
			
			
			// 7.2 mean
			gsl_vector_memcpy(mean_model[i], mean_prior[i]);
			gsl_vector_scale(mean_model[i], betai);
			
			gsl_vector_memcpy(temp, meank[i]);
			gsl_vector_scale(temp, nki);
			gsl_vector_add(mean_model[i], temp);
			gsl_vector_scale(mean_model[i], 1.0 / gsl_vector_get(beta_model, i));
			
			//7.3 wishart
			if(nki == 0.0) {
				gsl_matrix_memcpy(wish_model[i], wish_prior[i]);
			} else {
				
				gsl_matrix_memcpy(decomp, wish_prior[i]);
				gsl_linalg_LU_decomp(decomp, perm, &ind);
				gsl_linalg_LU_invert(decomp, perm, temp1);
				
				gsl_matrix_memcpy(temp2, sk[i]);
				gsl_matrix_scale(temp2, nki);
				gsl_matrix_add(temp1, temp2);
				
				double coeff = (betai * nki) / (betai + nki);
				gsl_vector_memcpy(temp, meank[i]);
				gsl_vector_sub(temp, mean_prior[i]);
				
				gsl_blas_dger(coeff, temp, temp, temp1);
				gsl_linalg_LU_decomp(temp1, perm, &ind);
				gsl_linalg_LU_invert(temp1, perm, wish_model[i]);
				
			}
		}
				
				
		// 4.1 : compute quad
		// new quad block : only 1 tempmat and meandat objects
		
		// new block bis : no nxn matrix any more...
		for(i=0; i<n; i++) {
			gsl_vector_view view1 = gsl_matrix_row(dat, i);
			for(j=0; j<k; j++) {
				gsl_vector_memcpy(meanvec, &view1.vector);
				gsl_vector_sub(meanvec, mean_model[j]);
				gsl_blas_dsymv(CblasUpper, gsl_vector_get(nu_model, j), wish_model[j], meanvec, 0.0, temp);
				gsl_blas_ddot(temp, meanvec, &val);
				val += d / gsl_vector_get(beta_model, j);
				gsl_matrix_set(quad, i, j, val);
			}
		}
		
		
		
		

	
		
		// 4.2 : compute kappalog
		for(i=0; i<d; i++) {
			gsl_vector_view view1 = gsl_matrix_column(gammamat, i);
			
			for(j=0; j<k; j++) {
				val = (gsl_vector_get(nu_model, j) - i) / 2.0;
				gsl_vector_set(&view1.vector, j, digamma(val));
			}
		}
		
		
		for(i=0; i<k; i++) {
			gsl_vector_view view1 = gsl_matrix_row(gammamat, i);
			tab = gsl_vector_ptr(&view1.vector, 0);
			val = sum(tab, d, view1.vector.stride);
			val += d * gsl_sf_log(2.0);
			
			gsl_matrix_memcpy(decomp, wish_model[i]);
			gsl_linalg_LU_decomp(decomp, perm, &ind);
			
			val += gsl_linalg_LU_lndet(decomp);
			gsl_vector_set(kappalog, i, val);
		}
		
		// 4.3 : compute pilog
		tab = gsl_vector_ptr(alpha_model, 0);
		double pisum = sum(tab, k, alpha_model->stride);
		
		pisum = digamma(pisum);
		for(i=0; i<k; i++) {
			val = digamma(gsl_vector_get(alpha_model, i)) - pisum;
			gsl_vector_set(pilog, i, val);
		}



		// 5 computing statistics
		
		// 5.1 : building pimatrix and kappamatrix
		
		for(i=0; i<n; i++) {
			gsl_vector_view view1 = gsl_matrix_row(pimatrix, i);
			gsl_vector_memcpy(&view1.vector, pilog);
			view1 = gsl_matrix_row(kappamatrix, i);
			gsl_vector_memcpy(&view1.vector, kappalog);
		}
		
		
		// 5.2 building resp
		// save resp for agitation calculation
		gsl_matrix_memcpy(oldresp, resp);
		
		gsl_matrix_set_zero(resp);
		gsl_matrix_add(resp, pimatrix);
	
		gsl_matrix_scale(kappamatrix, 0.5);
		gsl_matrix_scale(quad, 0.5);
		
		gsl_matrix_add(resp, kappamatrix);
		gsl_matrix_sub(resp, constmatrix);
		gsl_matrix_sub(resp, quad);
		
		
		for(i=0; i<n; i++) {
			gsl_vector_view view1 = gsl_matrix_row(resp, i);
			tab = gsl_vector_ptr(&view1.vector, 0);
			
			// instead, detect highest value in lnresp, and normalize by these value before computing exp
			// => no special case any more
			val = gsl_vector_max(&view1.vector);
			for(j=0; j<k; j++) {
				gsl_vector_set(&view1.vector, j, exp(gsl_vector_get(&view1.vector, j) - val));
			}
	
			val = sum(tab, k, view1.vector.stride);
			for(j=0; j<k; j++) {
				tab[j * view1.vector.stride] = tab[j * view1.vector.stride] / val;
			}
		}
		
		
		// 5.3 building statistics
				
		for(i=0; i<k; i++) {
			gsl_vector_view view1 = gsl_matrix_column(resp, i);
			tab = gsl_vector_ptr(&view1.vector, 0);
			val = sum(tab, n, view1.vector.stride);
			gsl_vector_set(nk, i, val);
			
			gsl_vector_set_zero(meank[i]);
			
			if(val != 0.0) {
				gsl_blas_dgemv(CblasTrans, 1.0/val, dat, &view1.vector, 0.0, meank[i]);
			}
			
			gsl_matrix_set_zero(sk[i]);
			
			if(val != 0.0) {
				for(j=0; j<n; j++) {
					gsl_vector_view view2 = gsl_matrix_row(dat, j);
					gsl_vector_memcpy(meanvec, &view2.vector);
					gsl_vector_sub(meanvec, meank[i]);
					gsl_blas_dsyr(CblasUpper, gsl_vector_get(&view1.vector, j), meanvec, sk[i]);
				}
				
				gsl_matrix_scale(sk[i], 1.0/gsl_vector_get(nk, i));
				uppercomplete(sk[i]);
			}
			
		}
		
		// 5.4 compute agitation
		for(i=0; i<k; i++) {
			agitation[i] = 0.0;
			if(gsl_vector_get(nk, i) > 0.01) {
				for(j=0; j<n; j++) {
					agitation[i] += fabs(gsl_matrix_get(resp, j, i) - gsl_matrix_get(oldresp, j, i));
				}
				agitation[i] = agitation[i] / gsl_vector_get(nk, i);
			}
		}
		

		// 6 compute bound value
		
		// 6.1 pxlog
		
		pxlog = 0.0;
		for(i=0; i<k; i++) {
			double tmp = 0.0;
			tmp += gsl_vector_get(kappalog, i);
			tmp -= d / gsl_vector_get(beta_model, i);
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, gsl_vector_get(nu_model, i), sk[i], wish_model[i], 0.0, prod);
			gsl_vector_view view1 = gsl_matrix_diagonal(prod);
			tab = gsl_vector_ptr(&view1.vector, 0);
			val = sum(tab, d, view1.vector.stride);
			tmp -= val;
			gsl_vector_memcpy(meanvec, meank[i]);
			gsl_vector_sub(meanvec, mean_model[i]);
			gsl_blas_dgemv(CblasNoTrans, gsl_vector_get(nu_model, i), wish_model[i], meanvec, 0.0, temp);
			gsl_vector_mul(temp, meanvec);
			tab = gsl_vector_ptr(temp, 0);
			val = sum(tab, d, temp->stride);
			tmp -= val;
			val = d * gsl_sf_log(2.0 * M_PI);
			tmp -= val;
			tmp = tmp * 0.5 * gsl_vector_get(nk, i);
			pxlog += tmp;
		}
			
		// 6.2 pzlog
		gsl_matrix_memcpy(calc, resp);
		gsl_matrix_mul_elements(calc, pimatrix);
		tab = gsl_matrix_ptr(calc, 0, 0);
		pzlog = sum(tab, n*k, 1);

		
		// 6.3 ppilog
		ppilog = 0.0;
		tab = gsl_vector_ptr(alpha_prior, 0);
		val = sum(tab, k, alpha_prior->stride);
		ppilog = gsl_sf_lngamma(val);
		for(i=0; i<k; i++) {
			ppilog -= gsl_sf_lngamma(gsl_vector_get(alpha_prior, i));
		}
		tab = gsl_vector_ptr(pilog, 0);
		val = sum(tab, k, pilog->stride);
		ppilog += (gsl_vector_get(alpha_prior, 0) - 1.0) * val;
		
		
		// 6.4 pmulog
		pmulog = 0.0;
		for(i=0; i<k; i++) {
			double tmp = 0.0;
			tmp += d * gsl_sf_log(gsl_vector_get(beta_prior, 0) / (2.0 * M_PI));
			tmp += gsl_vector_get(kappalog, i);
			tmp -= d * gsl_vector_get(beta_prior, 0) / gsl_vector_get(beta_model, i);
			
			// prod quadratique a priori non nul
			gsl_vector_memcpy(meanvec, mean_model[i]);
			gsl_vector_sub(meanvec, mean_prior[i]);
			gsl_blas_dgemv(CblasNoTrans, gsl_vector_get(beta_prior, 0) * gsl_vector_get(nu_model, i), wish_model[i], meanvec, 0.0, tempvec);
			gsl_vector_mul(tempvec, meanvec);
			tab = gsl_vector_ptr(tempvec, 0);
			tmp -= sum(tab, d, tempvec->stride);
			
			pmulog += 0.5 * tmp;
		}
		
		
		
		gsl_matrix_memcpy(decomp, wish_prior[0]);
		gsl_linalg_LU_decomp(decomp, perm, &ind);
			
		val = -(gsl_vector_get(nu_prior, 0) / 2.0) * gsl_linalg_LU_lndet(decomp);
		val -= gsl_vector_get(nu_prior, 0) * d * gsl_sf_log(2) / 2.0;
		val -= d * (d-1) * gsl_sf_log(M_PI) / 4.0;
		for(i=0; i<d; i++) {
			val -= gsl_sf_lngamma((gsl_vector_get(nu_prior, 0) - i) / 2.0);
		}
		
		pmulog += k * val;
		
		
		
		
		tab = gsl_vector_ptr(kappalog, 0);
		val = sum(tab, k, kappalog->stride);
		
		pmulog += (gsl_vector_get(nu_prior, 0) - d - 1.0) * val / 2.0;
		
		
		
		for(i=0; i<k; i++) {
			
			gsl_matrix_memcpy(decomp, wish_prior[i]);
			gsl_linalg_LU_decomp(decomp, perm, &ind);
			gsl_linalg_LU_invert(decomp, perm, inverse);
			
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverse, wish_model[i], 0.0, decomp);
			gsl_vector_view view1 = gsl_matrix_diagonal(decomp);
			tab = gsl_vector_ptr(&view1.vector, 0);
			
			pmulog -= 0.5 * gsl_vector_get(nu_model, i) * sum(tab, d, view1.vector.stride);
			
		}
		
		// 6.5 qzlog
		qzlog = 0.0;
		tab = gsl_matrix_ptr(resp, 0, 0);
		for(i=0; i<=(n-1)*(k-1); i++) {
			if(tab[i] != 0) {
				qzlog += tab[i] * gsl_sf_log(tab[i]);
			}
		}


		
		// 6.6 qpilog
		qpilog = 0.0;
		for(i=0; i<k; i++) {
			qpilog += (gsl_vector_get(alpha_model, i) - 1.0) * gsl_vector_get(pilog, i);
		}
		
		
		tab = gsl_vector_ptr(alpha_model, 0);
		val = sum(tab, k, alpha_model->stride);
		qpilog += gsl_sf_lngamma(val);
		for(i=0; i<k; i++) {
			qpilog -= gsl_sf_lngamma(gsl_vector_get(alpha_model, i));
		}
		
		// 6.7 qmulog
		qmulog = 0.0;
	
		
		for(i=0; i<k; i++) {
			double L = 0.0;
			double blog;
			for(j=0; j<d; j++) {
				L += gsl_sf_psi((gsl_vector_get(nu_model, i) - j) / 2.0);
			}
			gsl_matrix_memcpy(decomp, wish_model[i]);
			gsl_linalg_LU_decomp(decomp, perm, &ind);
			L += gsl_linalg_LU_lndet(decomp) + d * gsl_sf_log(2.0);
			
			blog = -(gsl_vector_get(nu_model, i) / 2.0) * gsl_linalg_LU_lndet(decomp);
			blog -= gsl_vector_get(nu_model, i) * d * gsl_sf_log(2) / 2.0;
			blog -= d * (d-1) * gsl_sf_log(M_PI) / 4.0;
			for(j=0; j<d; j++) {
				blog -= gsl_sf_lngamma((gsl_vector_get(nu_model, i) - j) / 2.0);
			}
			
			val = -((gsl_vector_get(nu_model, i) - d - 1.0) / 2.0) * L + (gsl_vector_get(nu_model, i) * d) / 2 - blog;
			
			qmulog += 0.5 * gsl_vector_get(kappalog, i) + d * gsl_sf_log(gsl_vector_get(beta_model, i) / (2.0 * M_PI)) / 2.0 - d / 2.0 - val;
			
		}
		
		// 8 : check convergence
		
		newbound = pxlog + pzlog + ppilog + pmulog - qzlog - qpilog - qmulog;
		
		count++;
		Rprintf("it %d bound is %f\n", count, newbound);
		
		// store nk and agitation
		if(count == 2) {
			nks = Calloc(1, double *);
			agitations = Calloc(1, double *);
		} else {
			nks = Realloc(nks, count-1, double *);
			agitations = Realloc(agitations, count-1, double *);
		}
		nks[count-2] = Calloc(k, double);
		agitations[count-2] = Calloc(k, double);
		for(i=0; i<k; i++) {
			nks[count-2][i] = gsl_vector_get(nk, i);
			agitations[count-2][i] = agitation[i];
		}



		// 2 kinds of convergence control : threshold, or if specified a number of iterations.
		if(maxit == R_NilValue) {
			if(newbound - bound < threshold) {
				converge = 1;
			}
		} else {
			if(count >= INTEGER(coerceVector(maxit, INTSXP))[0]) {
				converge = 1;
			}
		}
		
		bound = newbound;
	}

	
	// 9 : build return objects
	SEXP ret, mod, al, be, nu, me, wi, re;
	SEXP mecur, wicur, dims1, dims2;
	SEXP Rnks, Ragitations;
	
	PROTECT(ret = allocVector(VECSXP, 4));
	PROTECT(mod = allocVector(VECSXP, 6));
	PROTECT(al = allocVector(REALSXP, k));
	PROTECT(be = allocVector(REALSXP, k));
	PROTECT(nu = allocVector(REALSXP, k));
	PROTECT(me = allocVector(VECSXP, k));
	PROTECT(wi = allocVector(VECSXP, k));
	
	PROTECT(Rnks = allocVector(VECSXP, count-2));
	PROTECT(Ragitations = allocVector(VECSXP, count-2));
	
	double *altab = REAL(al);
	double *betab = REAL(be);
	double *nutab = REAL(nu);
	double *metab;
	double *witab;
	
	for(i=0; i<k; i++) {
		altab[i] = gsl_vector_get(alpha_model, i);
		betab[i] = gsl_vector_get(beta_model, i);
		nutab[i] = gsl_vector_get(nu_model, i);
		PROTECT(mecur = allocVector(REALSXP, d));
		metab = REAL(mecur);
		for(j=0; j<d; j++) {
			metab[j] = gsl_vector_get(mean_model[i], j);
		}
		SET_VECTOR_ELT(me, i, mecur);
		UNPROTECT(1);
		
		PROTECT(wicur = allocMatrix(REALSXP, d, d));
		witab = REAL(wicur);
		for(j=0; j<d; j++) {
			for(m=0; m<d; m++) {
				witab[j + m*d] = gsl_matrix_get(wish_model[i], j, m);
			}
		}
		SET_VECTOR_ELT(wi, i, wicur);
		UNPROTECT(1);
	}
	
	PROTECT(re = allocMatrix(REALSXP, n, k));
	tab = REAL(re);
	for(i=0; i<n; i++) {
		for(j=0; j<k; j++) {
			tab[i + n*j] = gsl_matrix_get(resp, i, j);
		}
	}

	for(i=0; i<count-2; i++) {
		SET_VECTOR_ELT(Rnks, i, allocVector(REALSXP, k));
		SET_VECTOR_ELT(Ragitations, i, allocVector(REALSXP, k));
		for(j=0; j<k; j++) {
			REAL(VECTOR_ELT(Rnks, i))[j] = nks[i][j];
			REAL(VECTOR_ELT(Ragitations, i))[j] = agitations[i][j];
		}
	}
	
	SET_VECTOR_ELT(mod, 0, al);
	SET_VECTOR_ELT(mod, 1, be);
	SET_VECTOR_ELT(mod, 2, nu);
	SET_VECTOR_ELT(mod, 3, me);
	SET_VECTOR_ELT(mod, 4, wi);
	SET_VECTOR_ELT(mod, 5, re);
	

	
	
	PROTECT(dims1 = allocVector(VECSXP, 4));
	PROTECT(dims2 = allocVector(VECSXP, 6));
	SET_VECTOR_ELT(dims1, 0, mkChar("model"));
	SET_VECTOR_ELT(dims1, 1, mkChar("data"));
	SET_VECTOR_ELT(dims1, 2, mkChar("nk"));
	SET_VECTOR_ELT(dims1, 3, mkChar("agitation"));
	SET_VECTOR_ELT(dims2, 0, mkChar("alpha"));
	SET_VECTOR_ELT(dims2, 1, mkChar("beta"));
	SET_VECTOR_ELT(dims2, 2, mkChar("nu"));
	SET_VECTOR_ELT(dims2, 3, mkChar("mean"));
	SET_VECTOR_ELT(dims2, 4, mkChar("wish"));
	SET_VECTOR_ELT(dims2, 5, mkChar("resp"));
	
	setAttrib(mod, R_NamesSymbol, dims2);
	
	SET_VECTOR_ELT(ret, 0, mod);
	SET_VECTOR_ELT(ret, 1, data);
	SET_VECTOR_ELT(ret, 2, Rnks);
	SET_VECTOR_ELT(ret, 3, Ragitations);
	
	setAttrib(ret, R_NamesSymbol, dims1);
	
			
		
		

	// 10 : free all alloc memory before return
	gsl_matrix_free(dat);
	gsl_matrix_free(quad);
	gsl_vector_free(pilog);
	gsl_vector_free(kappalog);
	gsl_matrix_free(resp);
	gsl_matrix_free(oldresp);
	gsl_vector_free(nk);
	for(i=0; i<k; i++) {
		gsl_vector_free(meank[i]);
	}
	free(meank);
	for(i=0; i<k; i++) {
		gsl_matrix_free(sk[i]);
	}
	free(sk);
	gsl_vector_free(alpha_prior);
	gsl_vector_free(beta_prior);
	gsl_vector_free(nu_prior);
	for(i=0; i<k; i++) {
		gsl_vector_free(mean_prior[i]);
		gsl_matrix_free(wish_prior[i]);
	}
	free(mean_prior);
	free(wish_prior);

	gsl_vector_free(alpha_model);
	gsl_vector_free(beta_model);
	gsl_vector_free(nu_model);
	for(i=0; i<k; i++) {
		gsl_vector_free(mean_model[i]);
		gsl_matrix_free(wish_model[i]);
	}
	free(mean_model);
	free(wish_model);
	
	// 11 : free all temporary variables
	gsl_matrix_free(c);
	gsl_matrix_free(betamatrix);
	gsl_vector_free(comput1);
	gsl_matrix_free(gammamat);
	gsl_matrix_free(decomp);
	gsl_matrix_free(pimatrix);
	gsl_matrix_free(kappamatrix);
	gsl_matrix_free(constmatrix);
	gsl_matrix_free(meanmatrix);
	gsl_vector_free(temp);
	gsl_matrix_free(temp1);
	gsl_matrix_free(temp2);
	gsl_matrix_free(interm);
	gsl_matrix_free(prod);
	gsl_vector_free(meanvec);
	gsl_matrix_free(calc);
	gsl_matrix_free(inverse);
	gsl_permutation_free(perm);
	gsl_vector_free(tempvec);
	gsl_vector_free(templine);
	
	for(i=0; i<count-2; i++) {
		Free(nks[i]);
		Free(agitations[i]);
	}
	Free(nks);
	Free(agitations);
	

	UNPROTECT(15);
	return(ret);

}
