// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include "utils.h"

void toLower(double *source, double *dest, int n) {
	for(int i=0; i<n; i++) {
		for(int j=0; j<=i; j++) {
			// !!! use n-1 instead of n (see reference in inst/lowerpack.pdf)
			dest[i+j*(1-j+2*(n-1))/2] = source[j*n+i];
		}
	}
}


void printSxpIntVector(SEXP sxp) {
	// assume double values
	int n = length(sxp);
	int *read = INTEGER(sxp);
	for(int i=0; i<n; i++) {
		Rprintf("%d ", read[i]);
	}
	Rprintf("\n");
}

void printSxpRealVector(SEXP sxp) {
	// assume double values
	int n = length(sxp);
	double *read = REAL(sxp);
	for(int i=0; i<n; i++) {
		Rprintf("%f ", read[i]);
	}
	Rprintf("\n");
}

void printSxpMatrix(SEXP sxp) {
	int n = INTEGER(getAttrib(sxp, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(sxp, R_DimSymbol))[1];
	double *read = REAL(sxp);
	for(int i=0; i<n; i++) {
		for(int j=0; j<d; j++) {
			Rprintf("%f ", read[i + n*j]);
		}
		Rprintf("\n");
	}
}


int imin(double *vec, int length, int incr) {
	int ind = 0;
	for(int i=0; i<length; i++) {
		if(vec[ind*incr] > vec[i*incr]) ind = i;
	}
	return(ind);
}

int imax(double *vec, int length, int incr) {
	int ind = 0;
	for(int i=0; i<length; i++) {
		if(vec[ind*incr] < vec[i*incr]) ind = i;
	}
	return(ind);
}


double vmin(double *vec, int length, int incr) {
	double val=R_PosInf;
	for(int i=0; i<length; i++) {
		if(val > vec[i*incr]) val = vec[i*incr];
	}
	return(val);
}

double vmax(double *vec, int length, int incr) {
	double val=R_NegInf;
	for(int i=0; i<length; i++) {
		if(val < vec[i*incr]) val = vec[i*incr];
	}
	return(val);
}

int allequal(int *vec1, int *vec2, int len) {
	int res=1;
	for(int i=0; i<len; i++) {
		if(vec1[i] != vec2[i]) res=0;
	}
	return(res);
}


int whichmax(gsl_vector *vect) {
	int ind;
	int n = vect->size;
	double max = GSL_NEGINF;
	for(int i=0; i<n; i++) {
		if(gsl_vector_get(vect, i) > max) {
			ind = i;
			max = gsl_vector_get(vect, i);
		}
	}
	return(ind);
}


inline SEXP getListElement(SEXP list, const char *str) { 
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol); 
	int i; 
	for (i = 0; i < length(list); i++) {
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) { 
			elmt = VECTOR_ELT(list, i); 
			break; 
		} 
	}
	return elmt; 
}

void vectorToSXP(SEXP *sxp, gsl_vector *vect) {
	double *read = REAL(*sxp);
	int size = vect->size;
	int i;
	
	for(i=0;i<size;i++) {
		read[i] = gsl_vector_get(vect, i);
	}
}

void intVectorToSXP(SEXP *sxp, gsl_vector *vect) {
	int *read = INTEGER(*sxp);
	int size = vect->size;
	int i;
	
	for(i=0;i<size;i++) {
		read[i] = (int)gsl_vector_get(vect, i);
	}

}

void matrixToSXP(SEXP *sxp, gsl_matrix *mat) {
	double *read = REAL(*sxp);
	int size1 = mat->size1;
	int size2 = mat->size2;
	int i,j;
	
	for(i=0;i<size1;i++) {
		for(j=0;j<size2;j++) {
			read[i + size1*j] = gsl_matrix_get(mat, i, j);
		}
	}
}

void SXPtoVector(gsl_vector *vect, SEXP sxp) {
	double *read = REAL(sxp);
	int size = length(sxp);
	int i;
	
	for(i=0; i<size; i++) {
		gsl_vector_set(vect, i, read[i]);
	}
	
}



void SXPtoMatrix(gsl_matrix *mat, SEXP sxp) {
	double *read = REAL(sxp);
	int n = INTEGER(getAttrib(sxp, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(sxp, R_DimSymbol))[1];
	int i,j;
	
	for(i=0; i<n; i++) {
		for(j=0; j<d; j++) {
			gsl_matrix_set(mat, i, j, read[i + n*j]);
		}
	}
}

double SXPvectorSum(SEXP sxp) {
	int i;
	double tot=0.0;
	int n = length(sxp);
	for(i=0; i<n; i++) {
		tot += REAL(sxp)[i];
	}
	return(tot);
	
}



double GSLvectorSum(gsl_vector *vect) {
	int i;
	int n = vect->size;
	int stride = vect->stride;
	double tot=0.0;
	for(i=0; i<n; i++) {
		//tot += gsl_vector_get(vect, i);
		tot += vect->data[i*stride];
	}
	return(tot);
}


// returns position sampled from the weight vector (from 0 to n-1)
SEXP multinomial(SEXP weights, SEXP nitems) {
	SEXP res;
	PROTECT(nitems= coerceVector(nitems, INTSXP));
	PROTECT(weights = coerceVector(weights, REALSXP));
	
	PROTECT(res = allocVector(INTSXP, INTEGER(nitems)[0]));
	
	double *read = REAL(weights);
	int i;
	int n = INTEGER(nitems)[0];
	
	GetRNGstate();
	for(i=0; i<n; i++) {
		double samp = runif(0.0, 1.0);
		double sum = read[0];
		int ind=0;
		while(sum < samp) {
			ind++;
			sum += read[ind];
		}
		INTEGER(res)[i] = ind;
	}
	PutRNGstate();
	UNPROTECT(3);
	return res;
}

SEXP mvngen(SEXP modmean, SEXP modcov, SEXP nitems) {
	double *read;
	int *iread;
	
	PROTECT(modmean = coerceVector(modmean, REALSXP));
	PROTECT(modcov = coerceVector(modcov, REALSXP));
	PROTECT(nitems = coerceVector(nitems, INTSXP));
	
	int n = INTEGER(nitems)[0];
	int d = length(modmean);
	int i,j;
	
	SEXP res;
	PROTECT(res = allocMatrix(REALSXP, n, d));
	// allocate analogous gsl matrix for calculations
	gsl_matrix *tempres = gsl_matrix_calloc(n,d);
	gsl_matrix *tempres2 = gsl_matrix_calloc(n,d);
	
	gsl_vector *mean = gsl_vector_calloc(d);
	gsl_matrix *cov = gsl_matrix_calloc(d,d);
	
	SXPtoVector(mean, modmean);
	SXPtoMatrix(cov, modcov);
	
	// fill result matrix with normal samples
	GetRNGstate();
	for(i=0; i<n; i++) {
		for(j=0; j<d; j++) {
			gsl_matrix_set(tempres, i, j, rnorm(0,1));
		}
	}
	PutRNGstate();
	
	// cholesky transform the covariance
	gsl_linalg_cholesky_decomp(cov);
	// lower to upper
	gsl_matrix_transpose(cov);
	// zero lower triangle
	for(i=1; i<d; i++) {
		for(j=0; j<i; j++) {
			gsl_matrix_set(cov, i, j, 0.0);
		}
	}
	
	// multiplication
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempres, cov, 0.0, tempres2);
	
	// ajout mean à chaque ligne
	gsl_vector_view view1;
	
	for(i=0; i<n; i++) {
		view1 = gsl_matrix_row(tempres2, i);
		gsl_vector_add(&view1.vector, mean);
	}
	
	matrixToSXP(&res, tempres2);
	
	gsl_matrix_free(tempres);
	gsl_matrix_free(tempres2);
	gsl_matrix_free(cov);
	gsl_vector_free(mean);
	
	UNPROTECT(4);
	return(res);
	
}


SEXP gmmgen(SEXP mod, SEXP nitem) {
	// input : w-mean-cov structured model mod, number of components k = length(w),
	// number of dimensions d=length(mean[1])
	// output : nxd matrix of sampled points
	double *read;
	int *iread;
	SEXP res, mn, cv;
	SEXP curcomp, curk;
	PROTECT(mod = coerceVector(mod, VECSXP));
	PROTECT(nitem = coerceVector(nitem, INTSXP));
	
	int n = INTEGER(nitem)[0];
	
	int i,j,l;
	int k,d=0;
	int cur;
	// w
	curcomp = coerceVector(getListElement(mod, "w"), REALSXP);
	k=length(curcomp);
	
	gsl_vector *w = gsl_vector_calloc(k);
	gsl_vector **mean = Calloc(k, gsl_vector *);
	gsl_matrix **cov = Calloc(k, gsl_matrix *);

	read = REAL(curcomp);	
	for(i=0; i<k; i++) {
		gsl_vector_set(w, i, read[i]);
	}		

	// mean
	curcomp = coerceVector(getListElement(mod, "mean"), VECSXP);
	for(i=0; i<k; i++) {
		curk = coerceVector(VECTOR_ELT(curcomp, i), REALSXP);
		read = REAL(curk);
		
		if(d==0) d = length(curk);
		
		mean[i] = gsl_vector_calloc(d);
		for(j=0; j<d; j++) {
			gsl_vector_set(mean[i], j, read[j]);
		}
	}
	
	// cov
	curcomp = coerceVector(getListElement(mod, "cov"), VECSXP);
	for(i=0; i<k; i++) {
		curk = coerceVector(VECTOR_ELT(curcomp, i), REALSXP);
		read = REAL(curk);
		
		cov[i] = gsl_matrix_calloc(d,d);
		for(j=0; j<d; j++) {
			for(l=0; l<d; l++) {
				gsl_matrix_set(cov[i], j, l, read[j + d*l]);
			}
		}
	}
	
	
	// 2 steps sampling
	// before proceeding, allocate matrix
	PROTECT(res=allocMatrix(REALSXP, n, d));
	// allocate objects to be passed to mvngen
	PROTECT(mn = allocVector(REALSXP, d));
	PROTECT(cv = allocMatrix(REALSXP, d,d));

	SEXP workvect;
	SEXP workvect2;
	PROTECT(workvect = allocVector(REALSXP, k));
	SEXP numb;
	PROTECT(numb = allocVector(INTSXP, 1));

	read = REAL(res);
	SEXP inds;
	PROTECT(inds = allocVector(INTSXP, n));
	
	// optimisation : 
	// first sample multinomial vector,
	// then get indexes of concerned occurences in a double vector,
	// and finally sample make a loop from 1 to k.
	// use lists to finely handle multinomial case
	SEXP multivec;
	vectorToSXP(&workvect, w);
	PROTECT(multivec = coerceVector(multinomial(workvect, nitem), INTSXP));
	
	list **root = Calloc(k, list *);
	for(i=0; i<k; i++) {
		//root[i] = createList();
		root[i] = Calloc(1, list);
	}
	
	for(i=0; i<n; i++) {
		cur=INTEGER(multivec)[i];
		addItem(root[cur], i, 0.0);
	}
	
	//for(i=0; i<k; i++) {
	//	for(j=0; j<size(root[i]); j++) {
	//		Rprintf("%d ", getItem(root[i], j));
	//	}
	//	Rprintf("\n");
	//}
	
	int offset=0;
	double *read2;
	
	for(i=0; i<k; i++) {		
		vectorToSXP(&mn, mean[i]);
		matrixToSXP(&cv, cov[i]);
		
		if(size(root[i]) != 0) {
			INTEGER(numb)[0] = size(root[i]);
			workvect2 = coerceVector(mvngen(mn, cv, numb), REALSXP);
			read2 = REAL(workvect2);
		}
		
		
		for(j=0; j<size(root[i]); j++) {
			for(l=0; l<d; l++) {
				read[offset + j + l*n] = read2[j + l*size(root[i])];
			}
			INTEGER(inds)[offset+j] = i+1;
		}
		
		offset += size(root[i]);
			
	}
	
	
	
	gsl_vector_free(w);
	for(i=0; i<k; i++) {
		gsl_vector_free(mean[i]);
		gsl_matrix_free(cov[i]);
		dropList(root[i]);
		Free(root[i]);
	}
	Free(mean);
	Free(cov);
	Free(root);
	
	SEXP pool;
	PROTECT(pool = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(pool, 0, res);
	SET_VECTOR_ELT(pool, 1, inds);
	
	UNPROTECT(10);
	return(pool);
		
			
}

// multivariate normal density of a datamatrix
// refactored to support rescaling, and not use GSL
SEXP mvndensity(SEXP mean, SEXP cov, SEXP data, SEXP rescaled) {
	int n = INTEGER(getAttrib(data, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(data, R_DimSymbol))[1];	
	SEXP res;
	PROTECT(res=allocVector(REALSXP, n));

	double *c_res = REAL(res);
	double *c_mean = REAL(mean);
	double *c_cov = REAL(cov);
	double *c_ccopy;
	double *c_data = REAL(data);
	int c_rescaled = INTEGER(rescaled)[0];
	double *c_samp = calloc(d, sizeof(double));
	double c_lndet;
	double *c_inverse;
	int i_one = 1;
	int dp1 = d+1;
	int d2 = d*d;

	// detect type in "general", "diagonal" and "spherical"
	char *type;
	double refval;
	int found=0;
	for(int i=0; i<d; i++) {
		for(int j=0; j<i; j++) {
			if(c_cov[j*d+i] != 0.0) found = 1;
		}
	}

	if(found) {
		type = "general";
		c_ccopy = calloc(d*d, sizeof(double));
		F77_CALL(dcopy)(&d2, c_cov, &i_one, c_ccopy, &i_one);
		c_inverse = calloc(d*d, sizeof(double));
	} else {
		refval = c_cov[0];
		for(int i=1; i<d; i++) {
			if(c_cov[i*d+i] != refval) found = 1;
		}
		if(found) {
			type = "diagonal";
			c_ccopy = calloc(d, sizeof(double));
			F77_CALL(dcopy)(&d, c_cov, &dp1, c_ccopy, &i_one);
			c_inverse = calloc(d, sizeof(double));
		} else {
			type = "spherical";
			c_ccopy = calloc(1, sizeof(double));
			c_ccopy[0] = c_cov[0];
			c_inverse = calloc(1, sizeof(double));
		}
	}

	// perform pre-decompositions
	symdecomp(c_ccopy, c_inverse, &c_lndet, d, type);

	for(int i=0; i<n; i++) {
		F77_CALL(dcopy)(&d, c_data+i, &n, c_samp, &i_one);
		if(c_rescaled == TRUE) {
			c_res[i] = drmnorm(c_samp, c_mean, c_lndet, c_inverse, d, 0, type);
		} else {
			c_res[i] = dmnorm(c_samp, c_mean, c_lndet, c_inverse, d, 0, type);
		}
	}
	
	free(c_samp);
	free(c_inverse);
	UNPROTECT(1);
	return(res);
	
}


SEXP mvnradiusdensity(SEXP cov, SEXP radii) {
	
	int d = INTEGER(getAttrib(cov, R_DimSymbol))[0];
	int n = length(radii);
	int i;
	double val;
	gsl_permutation *perm = gsl_permutation_alloc(d);
	
	gsl_matrix *cv = gsl_matrix_calloc(d,d);
	gsl_vector *rad = gsl_vector_calloc(n);
	
	SXPtoMatrix(cv, cov);
	SXPtoVector(rad, radii);
	
	// return result in a vector
	gsl_vector *result = gsl_vector_calloc(n);
	
	// coeff and inverse matrix are the same in all cases
	gsl_matrix *decomp = gsl_matrix_calloc(d,d);
	gsl_matrix_memcpy(decomp, cv);
	gsl_linalg_LU_decomp(decomp, perm, &i);
	
	// compute log probs to prevent singularities
	double lndet = gsl_linalg_LU_lndet(decomp);
	gsl_matrix *inverse = gsl_matrix_alloc(d,d);
	gsl_linalg_LU_invert(decomp, perm, inverse);
	
	double constant = 0.0;
	constant += - (d / 2.0) * gsl_sf_log(2 * M_PI) - 0.5 * lndet;

	
	//gsl_vector *meanvec = gsl_vector_alloc(d);
	//gsl_vector *tempvec = gsl_vector_alloc(d);
	//gsl_vector_view view1;
	
	for(i=0; i<n; i++) {
		//view1 = gsl_matrix_row(dat, i);
		//gsl_vector_memcpy(meanvec, &view1.vector);
		//gsl_vector_sub(meanvec, mn);
		//
		//gsl_blas_dsymv(CblasUpper, -0.5, inverse, meanvec, 0.0, tempvec);
		//gsl_blas_ddot(meanvec, tempvec, &val);
		
		gsl_vector_set(result, i, exp(constant - 0.5 * gsl_vector_get(rad,i)));
	}
	
	SEXP res;
	PROTECT(res=allocVector(REALSXP, n));
	vectorToSXP(&res, result);
	
	gsl_permutation_free(perm);
	gsl_matrix_free(cv);
	gsl_vector_free(rad);
	gsl_vector_free(result);
	gsl_matrix_free(decomp);
	gsl_matrix_free(inverse);
	//gsl_vector_free(meanvec);
	//gsl_vector_free(tempvec);
	
	UNPROTECT(1);
	return(res);
	
}

void symdecomp(double *mat, double *inverse, double *lndet, int n, const char *type) {
	// from an input symmetric positive definite matrix (typically a covariance matrix)
	// compute summaries that accelerate dmnorm and likelihood computations

	// added optimization for diagonal matrices

	int i_one=1;
	int size;
	if(!strcmp(type, "general")) {
		size = n*n;
	} else if(!strcmp(type, "diagonal")) {
		size = n;
	} else {
		size = 1;
	}
	int info=0;
	
	F77_CALL(dcopy)(&size, mat, &i_one, inverse, &i_one);
	// compute Cholesky only in general case
	if(!strcmp(type, "general")) {
		F77_CALL(dpotrf)("L", &n, inverse, &n, &info);
	}
	
	lndet[0] = 0.0;
	for(int j=0; j<n; j++) {
		if(!strcmp(type, "general")) {
			lndet[0] = lndet[0] + log(inverse[j*n+j]);
		} else if(!strcmp(type, "diagonal")) {
			lndet[0] = lndet[0] + log(inverse[j]);
		} else {
			lndet[0] = lndet[0] + log(inverse[0]);
		}
	}

	if(!strcmp(type, "general")) {
		// if det computed from Cholesky
		lndet[0] = lndet[0] * 2.0;
	}

	if(!strcmp(type, "general")) {
		F77_CALL(dpotri)("L", &n, inverse, &n, &info);
	} else if(!strcmp(type, "diagonal")) {
		for(int j=0; j<n; j++) inverse[j] = 1.0 / inverse[j];
	} else {
		inverse[0] = 1.0 / inverse[0];
	}
	
	// make symmetric if type="general"
	// FIX: useless if only lower part used through calculations
	//if(!strcmp(type, "general")) {
	//	for(int i=0; i<n; i++) {
	//		for(int j=(i+1); j<n; j++) {
	//			inverse[j*n+i] = inverse[i*n+j];
	//		}
	//	}		
	//}

}


double dmnorm(double *samp, double *mean, double lndetcov, double *inverse, int d, int logd, const char *type) {
	// fast implementation using pre-computed lndet and inverse cov, to allow incremental computation
	// rather use cholesky decomp -> then just Ax and ddot suffices
	double *vec = calloc(d, sizeof(double));
	double *vec2 = calloc(d, sizeof(double));
	int i_one = 1;
	double alpha = -1.0;
	double beta = 0.0;
	F77_CALL(dcopy)(&d, samp, &i_one, vec, &i_one);
	F77_CALL(daxpy)(&d, &alpha, mean, &i_one, vec, &i_one);

	alpha = 1.0;
	double res=0.0;
	if(!strcmp(type, "general")) {
		F77_CALL(dsymv)("L", &d, &alpha, inverse, &d, vec, &i_one, &beta, vec2, &i_one);
		res = F77_CALL(ddot)(&d, vec, &i_one, vec2, &i_one);
	} else if(!strcmp(type, "diagonal")) {
		for(int i=0; i<d; i++) res += inverse[i] * pow(vec[i], 2.0);
	} else {
		for(int i=0; i<d; i++) res += inverse[0] * pow(vec[i], 2.0);
	}

	// compute in log to prevent singularities on log(det)
	res = -d*log(2.0*M_PI)/2.0 -0.5*lndetcov -0.5*res;
	if(!logd) {
		res = exp(res);
	}
	
	free(vec);
	free(vec2);
	return(res);
}

// rescaled unnormalized dnorm, so that BIC criterion is insensitive to dimensionality
// -> no 2PI const
// -> d power on detcov to limit tendency of det to explode exponentially
// -> -d to the squared norm in the exp() term, as is shown to be the expectation of squared radius
double drmnorm(double *samp, double *mean, double lndetcov, double *inverse, int d, int logd, const char *type) {
	double *vec = calloc(d, sizeof(double));
	double *vec2 = calloc(d, sizeof(double));
	int i_one = 1;
	double alpha = -1.0;
	double beta = 0.0;
	F77_CALL(dcopy)(&d, samp, &i_one, vec, &i_one);
	F77_CALL(daxpy)(&d, &alpha, mean, &i_one, vec, &i_one);

	alpha = 1.0;
	double res=0.0;
	if(!strcmp(type, "general")) {
		F77_CALL(dsymv)("L", &d, &alpha, inverse, &d, vec, &i_one, &beta, vec2, &i_one);
		res = F77_CALL(ddot)(&d, vec, &i_one, vec2, &i_one);
	} else if(!strcmp(type, "diagonal")) {
		for(int i=0; i<d; i++) res += inverse[i] * pow(vec[i], 2.0);
	} else {
		for(int i=0; i<d; i++) res += inverse[0] * pow(vec[i], 2.0);
	}

	// compute in log to prevent singularities on log(det)
	res = -(1.0/(2.0*d))*lndetcov - 0.5*(res-d);
	if(!logd) {
		res = exp(res);
	}
	
	free(vec);
	free(vec2);
	return(res);	
}



// gmm density of a data matrix
SEXP gmmdensity(SEXP mod, SEXP data) {
	
	
	SEXP curvect;
	SEXP totvect;
	SEXP r_false = PROTECT(allocVector(INTSXP, 1));
	INTEGER(r_false)[0] = 0;

	int n = INTEGER(getAttrib(data, R_DimSymbol))[0];
	int k = length(getListElement(mod, "w"));
	int i;

	
	gsl_vector *curgslvect = gsl_vector_calloc(n);
	gsl_vector *totgslvect = gsl_vector_calloc(n);
	
	for(i=0; i<k; i++) {
		curvect = mvndensity(VECTOR_ELT(getListElement(mod, "mean"), i), VECTOR_ELT(getListElement(mod, "cov"), i), data, r_false);
		SXPtoVector(curgslvect, curvect);
		gsl_vector_scale(curgslvect, REAL(getListElement(mod, "w"))[i] );
		gsl_vector_add(totgslvect, curgslvect);
	}
		
	PROTECT(totvect = allocVector(REALSXP, n));
	vectorToSXP(&totvect, totgslvect);
	
	gsl_vector_free(curgslvect);
	gsl_vector_free(totgslvect);
	
	UNPROTECT(2);
	return(totvect);
	
}

SEXP klmc(SEXP mod1, SEXP mod2, SEXP nit) {
	
	
	// sample from gmm mod1
	SEXP sample;
	PROTECT(sample = VECTOR_ELT(gmmgen(coerceVector(mod1, VECSXP), coerceVector(nit, INTSXP)), 0));

	int n = INTEGER(coerceVector(nit, INTSXP))[0];
	int i,j;
	double *read;
	
	
	
	// obtain density vectors 
	SEXP density1, density2;
	PROTECT(density1 = gmmdensity(coerceVector(mod1, VECSXP), sample));
	PROTECT(density2 = gmmdensity(coerceVector(mod2, VECSXP), sample));
	
	// convert to gsl vectors for calculation
	gsl_vector *gslden1 = gsl_vector_calloc(n);
	gsl_vector *gslden2 = gsl_vector_calloc(n);
	
	SXPtoVector(gslden1, density1);
	SXPtoVector(gslden2, density2);
	
	for(i=0; i<n; i++) {
		gsl_vector_set(gslden1, i, (1.0 / n) * gsl_sf_log(gsl_vector_get(gslden1, i) / gsl_vector_get(gslden2, i)));
	}
	
	double res = GSLvectorSum(gslden1);
	
	SEXP sxpRes;
	PROTECT(sxpRes = allocVector(REALSXP, 1));
	REAL(sxpRes)[0] = res;
	
	gsl_vector_free(gslden1);
	gsl_vector_free(gslden2);
	
	UNPROTECT(4);
	return(sxpRes);
	
}


SEXP jsmc(SEXP mod1, SEXP mod2, SEXP nit) {
	// Jensen shannon divergence between mod1 et mod2
	
	// build average model
	// size k1 + k2, with all coeffs half from original value
	
	PROTECT(mod1 = coerceVector(mod1, VECSXP));
	PROTECT(mod2 = coerceVector(mod2, VECSXP));
	PROTECT(nit = coerceVector(nit, INTSXP));
	
	int k1 = length(coerceVector(getListElement(mod1, "w"), REALSXP));
	int k2 = length(coerceVector(getListElement(mod2, "w"), REALSXP));
	int i;
	double *read;
	
	// build global w vector
	gsl_vector *globw = gsl_vector_alloc(k1+k2);
	read = REAL(coerceVector(getListElement(mod1, "w"), REALSXP));
	
	for(i=0; i<k1; i++) {
		gsl_vector_set(globw, i, read[i]);
	}
	
	read = REAL(coerceVector(getListElement(mod2, "w"), REALSXP));
	
	for(i=0; i<k2; i++) {
		gsl_vector_set(globw, k1+i, read[i]);
	}
	
	gsl_vector_scale(globw, 1.0/2.0);
	
	SEXP avmod, avw, avmean, avcov;
	
	PROTECT(avw = allocVector(REALSXP, k1+k2));
	vectorToSXP(&avw, globw);
	
	// build SXP target structure
	PROTECT(avmod = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(avmod, 0, avw);
	
	PROTECT(avmean = allocVector(VECSXP, k1+k2));
	PROTECT(avcov = allocVector(VECSXP, k1+k2));
	
	for(i=0; i<k1; i++) {
		SET_VECTOR_ELT(avmean, i, VECTOR_ELT(coerceVector(getListElement(mod1, "mean"), VECSXP), i));
		SET_VECTOR_ELT(avcov, i, VECTOR_ELT(coerceVector(getListElement(mod1, "cov"), VECSXP), i));
	}
	
	for(i=0; i<k2; i++) {
		SET_VECTOR_ELT(avmean, k1+i, VECTOR_ELT(coerceVector(getListElement(mod2, "mean"), VECSXP), i));
		SET_VECTOR_ELT(avcov, k1+i, VECTOR_ELT(coerceVector(getListElement(mod2, "cov"), VECSXP), i));
	}
	
	SET_VECTOR_ELT(avmod, 1, avmean);
	SET_VECTOR_ELT(avmod, 2, avcov);
	
	SEXP dims;
	PROTECT(dims = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(dims, 0, mkChar("w"));
	SET_VECTOR_ELT(dims, 1, mkChar("mean"));
	SET_VECTOR_ELT(dims, 2, mkChar("cov"));
	
	setAttrib(avmod, R_NamesSymbol, dims);
	
	double res = 0.5 * REAL(klmc(mod1, avmod, nit))[0] + 0.5 * REAL(klmc(mod2, avmod, nit))[0];
	
	SEXP sxpRes;
	PROTECT(sxpRes = allocVector(REALSXP, 1));
	REAL(sxpRes)[0] = res;
	
	gsl_vector_free(globw);
	
	UNPROTECT(9);
	return(sxpRes);
	
	
}

SEXP klut(SEXP mod1, SEXP mod2) {
	// computes klut(mod1 || mod2)
	// for 2 models structured (w, mean, cov)
	// assumes coherent dimensionalities
	PROTECT(mod1 = coerceVector(mod1, VECSXP));
	PROTECT(mod2 = coerceVector(mod2, VECSXP));

	int k1 = length(getListElement(mod1, "w"));
	int k2 = length(getListElement(mod2, "w"));
	int d = length(VECTOR_ELT(getListElement(mod1, "mean"), 0));
	int i,j;
	// create data structure for the sample
	// NB easier to have a sample per k1 components.
	// rather use Calloc
	gsl_matrix **samp = Calloc(k1 , gsl_matrix *);
	for(i=0; i<k1; i++) {
		samp[i] = gsl_matrix_alloc(2 * d + 1, d);
	}
	//gsl_matrix *samp = gsl_matrix_alloc(k1 * (2 * d + 1), d);
	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(4 * d);
	gsl_matrix *a = gsl_matrix_alloc(d,d);
	gsl_matrix *vecs = gsl_matrix_alloc(d,d);
	gsl_vector *vals = gsl_vector_alloc(d);
	gsl_vector *mean = gsl_vector_alloc(d);
	gsl_vector_view view1;
	gsl_vector_view view2;
	SEXP cursamp;
	PROTECT_INDEX ipx[2];
	PROTECT(cursamp = allocMatrix(REALSXP, 2*d+1, d));
	SEXP densities1;
	PROTECT_WITH_INDEX(densities1 = allocVector(REALSXP, 2*d+1), &ipx[0]);
	SEXP densities2;
	PROTECT_WITH_INDEX(densities2 = allocVector(REALSXP, 2*d+1), &ipx[1]);
	
		
	// perform eigen value decompositions, and fill sample matrix with its results
	for(i=0; i<k1; i++) {
		SXPtoMatrix(a, VECTOR_ELT(getListElement(mod1, "cov"), i));
		SXPtoVector(mean, VECTOR_ELT(getListElement(mod1, "mean"), i));
		gsl_eigen_symmv(a, vals, vecs, workspace);
		for(int j=0; j<d; j++) {
			view1 = gsl_matrix_row(samp[i], 2*j);
			view2 = gsl_matrix_column(vecs, j);
			gsl_vector_scale(&view2.vector, sqrt(gsl_vector_get(vals, j)));
			gsl_vector_memcpy(&view1.vector, mean);
			gsl_vector_add(&view1.vector, &view2.vector);
			
			view1 = gsl_matrix_row(samp[i], 2*j+1);
			gsl_vector_memcpy(&view1.vector, mean);
			gsl_vector_sub(&view1.vector, &view2.vector);
		}
		view1 = gsl_matrix_row(samp[i], 2*d);
		gsl_vector_memcpy(&view1.vector, mean);
	}
		
	// calculate KL approximation
	double klapprox = 0.0;
	for(i=0; i<k1; i++) {
		matrixToSXP(&cursamp, samp[i]);
		REPROTECT(densities1 = gmmdensity(mod1, cursamp), ipx[0]);
		REPROTECT(densities2 = gmmdensity(mod2, cursamp), ipx[1]);
		
		
		for(j=0; j<(2*d+1); j++) {
			// manage potential case where numerator or denominator are 0
			double num, den;
			num = REAL(densities1)[j];
			den = REAL(densities2)[j];
			
			
			if(num == 0.0) {
				num = pow(10, -10);
			}
			if(den == 0.0) {
				den = pow(10, -10);
			}
			klapprox += (1.0 / (2.0 * d + 1.0)) * (REAL(getListElement(mod1, "w"))[i]) * gsl_sf_log(num/den);
		}
	}
	
	SEXP result;
	PROTECT(result = allocVector(REALSXP, 1));
	// no reason for this value to be negative.
	if(klapprox<0.0) {
		REAL(result)[0] = 0.0;
	} else {
		REAL(result)[0] = klapprox;
	}
	
	// deallocate everything
	for(i=0; i<k1; i++) {
		gsl_matrix_free(samp[i]);
	}
	Free(samp);
	gsl_eigen_symmv_free(workspace);
	gsl_matrix_free(a);
	gsl_matrix_free(vecs);
	gsl_vector_free(vals);
	gsl_vector_free(mean);
	UNPROTECT(6);
	return(result);
		
}


SEXP jsut(SEXP mod1, SEXP mod2) {
	// jensen shannon with unscented transform
	PROTECT(mod1 = coerceVector(mod1, VECSXP));
	PROTECT(mod2 = coerceVector(mod2, VECSXP));
	
	int k1 = length(coerceVector(getListElement(mod1, "w"), REALSXP));
	int k2 = length(coerceVector(getListElement(mod2, "w"), REALSXP));
	int i;
	double *read;
	
	// build global w vector
	gsl_vector *globw = gsl_vector_alloc(k1+k2);
	read = REAL(coerceVector(getListElement(mod1, "w"), REALSXP));


	for(i=0; i<k1; i++) {
		gsl_vector_set(globw, i, read[i]);
	}
	
	read = REAL(coerceVector(getListElement(mod2, "w"), REALSXP));
	
	for(i=0; i<k2; i++) {
		gsl_vector_set(globw, k1+i, read[i]);
	}
	
	gsl_vector_scale(globw, 1.0/2.0);
	
	SEXP avmod, avw, avmean, avcov;
	
	PROTECT(avw = allocVector(REALSXP, k1+k2));
	vectorToSXP(&avw, globw);
	
	// build SXP target structure
	PROTECT(avmod = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(avmod, 0, avw);
	
	PROTECT(avmean = allocVector(VECSXP, k1+k2));
	PROTECT(avcov = allocVector(VECSXP, k1+k2));
	
	for(i=0; i<k1; i++) {
		SET_VECTOR_ELT(avmean, i, VECTOR_ELT(coerceVector(getListElement(mod1, "mean"), VECSXP), i));
		SET_VECTOR_ELT(avcov, i, VECTOR_ELT(coerceVector(getListElement(mod1, "cov"), VECSXP), i));
	}
	
	for(i=0; i<k2; i++) {
		SET_VECTOR_ELT(avmean, k1+i, VECTOR_ELT(coerceVector(getListElement(mod2, "mean"), VECSXP), i));
		SET_VECTOR_ELT(avcov, k1+i, VECTOR_ELT(coerceVector(getListElement(mod2, "cov"), VECSXP), i));
	}
	
	SET_VECTOR_ELT(avmod, 1, avmean);
	SET_VECTOR_ELT(avmod, 2, avcov);
	

	SEXP dims;
	PROTECT(dims = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(dims, 0, mkChar("w"));
	SET_VECTOR_ELT(dims, 1, mkChar("mean"));
	SET_VECTOR_ELT(dims, 2, mkChar("cov"));
	
	setAttrib(avmod, R_NamesSymbol, dims);
	
	//double res = 0.5 * REAL(klut(mod1, avmod))[0] + 0.5 * REAL(klut(mod2, avmod))[0];
	// protect klut results	
	SEXP klmod1, klmod2;
	PROTECT(klmod1 = klut(mod1, avmod));
	PROTECT(klmod2 = klut(mod2, avmod));
	double res = 0.5 * REAL(klmod1)[0] + 0.5 * REAL(klmod2)[0];


	SEXP sxpRes;
	PROTECT(sxpRes = allocVector(REALSXP, 1));
	REAL(sxpRes)[0] = res;
	
	gsl_vector_free(globw);
	
	UNPROTECT(10);
	return(sxpRes);
}



	
		

//list *createList() {
//	list *theList = Calloc(1, list);
//}

void addItem(list *theList, int intVal, double doubleVal) {
	if(theList->first == 0) {
		theList->first = Calloc(1, item);
		theList->first->intVal = intVal;
		theList->first->doubleVal = doubleVal;
		theList->last = theList->first;
	} else {
		theList->last->next = Calloc(1, item);
		theList->last->next->intVal = intVal;
		theList->last->next->doubleVal = doubleVal;
		theList->last = theList->last->next;
	}
}

void addList(metalist *metalst, int key) {
	metaitem *it = Calloc(1, metaitem);
	it->lst = Calloc(1, list);
	it->key = key;
	if(metalst->first == 0) {
		metalst->first = it;
		metalst->last = metalst->first;
	} else {
		metalst->last->next = it;
		metalst->last = metalst->last->next;
	}
}

int metaSize(metalist *metalst) {
	int sz=0;
	metaitem *ptr = metalst->first;
	while(ptr != 0) {
		ptr = ptr->next;
		sz++;
	}
	
	return sz;
}

list *getList(metalist *metalst, int index) {
	int i=0;
	metaitem *ptr = metalst->first;
	while(i<index) {
		if(ptr == 0) {
			return NULL;
		} else {
			ptr = ptr->next;
		}
		i++;
	}
	
	// guaranteed here to point to an actual element
	return ptr->lst;
}

int getIntMetaItem(metalist *meta, int index) {
	int i=0;
	metaitem *ptr = meta->first;
	while(i<index) {
		if(ptr == 0) {
			return 0;
		} else {
			ptr = ptr->next;
		}
		i++;
	}
	
	return ptr->key;
}
	
	
	

int indexOfInt(list *lst, int key) {
	int res=-1;
	int i;
	for(i=0; i<size(lst); i++) {
		if(getIntItem(lst, i) == key) {
			res=i;
		}
	}
	
	return res;
	
}


int indexOfKey(metalist *metalst, int key) {
	int res=-1;
	int i;
	metaitem *it = metalst->first;
	for(i=0; i<metaSize(metalst); i++) {
		if(it->key == key) {
			res = i;
		}
		it = it->next;
	}
	
	return res;
}

SEXP listToSXP(list *lst) {
	int len = size(lst);
	int i;
	
	SEXP res;
	PROTECT(res = allocVector(INTSXP, len));
	
	for(i=0; i<len; i++) {
		INTEGER(res)[i] = getIntItem(lst, i);
	}
	
	UNPROTECT(1);
	return res;
}
	

void insertItem(list *theList, int intVal, double doubleVal, int pos) {
	// inserts values after pos item
	// if pos > size - 1, do as if pos = size - 1
	// at least 1 element should be in the list
	int len = size(theList);
	if(pos > len - 1) {
		pos = len - 1;
	}
	
	int i;
	item *cur = theList->first;
	item *temp;
	for(i=0; i<pos; i++) {
		cur = cur->next;
	}
	
	temp = Calloc(1, item);
	temp->intVal = intVal;
	temp->doubleVal = doubleVal;
	temp->next = cur->next;
	cur->next = temp;
}
	

int size(list *theList) {
	int l=0;
	item *ptr = theList->first;
	while(ptr !=0) {
		l++;
		ptr = ptr->next;
	}
	return(l);
}

int getIntItem(list *theList, int index) {
	int l=0;
	item *ptr = theList->first;
	while(l<index) {
		l++;
		ptr = ptr->next;
	}
	return(ptr->intVal);
}

double getDoubleItem(list *theList, int index) {
	int l=0;
	item *ptr = theList->first;
	while(l<index) {
		l++;
		ptr = ptr->next;
	}
	return(ptr->doubleVal);
}

void dropItem(list *theList, int index) {
	// at least 1 item should be in the list
	int i;
	item *cur = theList->first;
	item *temp;
	
	for(i=0; i<index; i++) {
		temp = cur;
		cur = cur->next;
	}
	
	if(index == 0) {
		theList->first = cur->next;
		Free(cur);
	} else {
		temp->next = cur->next;
		Free(cur);
	}
}
	

void dropList(list *theList) {
	item *cur;
	item *ptr = theList->first;
	while(ptr != 0) {
		cur = ptr->next;
		Free(ptr);
		ptr = cur;
	}
}

void dropMetalist(metalist *meta) {
	metaitem *ptr = meta->first;
	metaitem *cur;
	
	while(ptr != 0) {
		Free(ptr->lst);
		cur = ptr->next;
		Free(ptr);
		ptr = cur;
	}
}


void setDoubleItem(list *lst, int index, double value) {
	item *ptr = lst->first;
	int i=0;
	while(i<index) {
		i++;
		ptr = ptr->next;
	}
	
	ptr->doubleVal = value;
}
	
	

void printList(list *lst) {
	item *ptr = lst->first;
	while(ptr != NULL) {
		Rprintf("%d %f\n", ptr->intVal, ptr->doubleVal);
		//printf("%d %f\n", ptr->intVal, ptr->doubleVal);
		ptr = ptr->next;
	}
}

		
	
void listToIntVector(gsl_vector *v, list *l) {
	for(int i=0; i<size(l); i++) {
		gsl_vector_set(v, i, getIntItem(l, i));
	}	
}

void listToDoubleVector(gsl_vector *v, list *l) {
	for(int i=0; i<size(l); i++) {
		gsl_vector_set(v, i, getDoubleItem(l, i));
	}	
}
	

	
SEXP extractSimpleModel(SEXP model, SEXP labels) {
	// extract a w, mean, cov model from a raw varbayes model
	// model : the model to process (list with model and data attributes
	// labels : boolean indicating wether to build a labels vector from data
	
	
	// test pour expr logique => OK
	//int toto = INTEGER(coerceVector(labels, LGLSXP))[0];
	//Rprintf("%d\n", toto);
	
	int i, j, n, k, d, kbis, hasData, ind, cur;
	double cumul;
	hasData = INTEGER(coerceVector(labels, LGLSXP))[0];
	
	if(hasData) {
		// in this case model has a data component structured as a matrix, we get its n
		// d will be taken from the model itself
		n = INTEGER(getAttrib(getListElement(model, "data"), R_DimSymbol))[0];
	}
	
	// for now, labels info is not used further
	
	SEXP alpha = getListElement(getListElement(model, "model"), "alpha");
	
	k = length(alpha);
	d = INTEGER(getAttrib(VECTOR_ELT(getListElement(getListElement(model, "model"), "wish"), 0), R_DimSymbol))[0];
	cumul = SXPvectorSum(alpha);
	
	//list *selection = createList();
	list *selection = Calloc(1, list);	

	for(i=0; i<k; i++) {
		// extract component only if more than 2 points supported
		if((REAL(alpha)[i] - 0.01) > 2.0) {
			addItem(selection, i, 0.0);
		}
	}
	
	Rprintf("selected size : %d\n", size(selection));
	kbis = size(selection);

	// create data structures for w and cov as calculations are needed there.
	// for mean we directly declare a list that will be instantiated with coerceVector from input structure
	gsl_vector *gslW = gsl_vector_alloc(kbis);
	gsl_matrix **gslCov = Calloc(kbis, gsl_matrix *);
	for(i=0; i<kbis; i++) {
		gslCov[i] = gsl_matrix_alloc(d,d);
	}
	gsl_matrix *decomp = gsl_matrix_alloc(d,d);
	gsl_permutation *perm = gsl_permutation_alloc(d);
	SEXP mean;
	PROTECT(mean = allocVector(VECSXP, kbis));
	
	gsl_matrix **inverses = Calloc(kbis, gsl_matrix *);
	for(i=0; i<kbis; i++) {
		inverses[i] = gsl_matrix_alloc(d,d);
	}
	double lndets[kbis];
	
	for(i=0; i<kbis; i++) {
		cur = getIntItem(selection, i);
		gsl_vector_set(gslW, i, REAL(alpha)[cur]);
		SET_VECTOR_ELT(mean, i, coerceVector(VECTOR_ELT(getListElement(getListElement(model, "model"), "mean"), cur), REALSXP));
		SXPtoMatrix(gslCov[i], VECTOR_ELT(getListElement(getListElement(model, "model"), "wish"), cur));
		gsl_matrix_memcpy(decomp, gslCov[i]);
		gsl_matrix_scale(decomp, REAL(getListElement(getListElement(model, "model"), "nu"))[cur]);
		gsl_matrix_memcpy(inverses[i], decomp);
		gsl_linalg_LU_decomp(decomp, perm, &ind);
		gsl_linalg_LU_invert(decomp, perm, gslCov[i]);
		lndets[i] = gsl_linalg_LU_lndet(decomp);
	}
	
	// scale w and create final data structures
	cumul = GSLvectorSum(gslW);
	gsl_vector_scale(gslW, 1.0 / cumul);
	
	// at this step we have the model : we compute arg max component for each line of data
	SEXP labs;

	if(hasData) {
		double curdensity;
		double constant = - d * gsl_sf_log(2.0 * M_PI) / 2.0;
		double theMax;
		PROTECT(labs=allocVector(INTSXP, n));
		gsl_vector_view view1;
		gsl_vector_view view2;
		gsl_vector *temp1 = gsl_vector_alloc(d);
		gsl_vector *temp2 = gsl_vector_alloc(d);
		gsl_matrix_view matview1;
		matview1 = gsl_matrix_view_array(REAL(getListElement(model, "data")), d, n);
		
		for(i=0; i<n; i++) {
			theMax = GSL_NEGINF;
			for(j=0; j<kbis; j++) {
				view1 = gsl_matrix_column(&matview1.matrix, i);
				gsl_vector_memcpy(temp1, &view1.vector);
				view2 = gsl_vector_view_array(REAL(VECTOR_ELT(mean, j)), d);
				gsl_blas_daxpy(-1.0, &view2.vector, temp1);
				gsl_blas_dsymv(CblasUpper, -0.5, inverses[j], temp1, 0.0, temp2);
				gsl_blas_ddot(temp1, temp2, &curdensity);
				curdensity += constant - 0.5 * lndets[j];
				if(curdensity > theMax) {
					INTEGER(labs)[i] = j+1;
					theMax = curdensity;
				}
			}
		}
		
		gsl_vector_free(temp1);
		gsl_vector_free(temp2);
		
	}	
			
	
	
	SEXP w, cov, res;
	PROTECT(w = allocVector(REALSXP, kbis));
	PROTECT(cov = allocVector(VECSXP, kbis));
	if(hasData) {
		PROTECT(res = allocVector(VECSXP, 4));
	} else {
		PROTECT(res = allocVector(VECSXP, 3));
	}
	
	SEXP dims;
	if(hasData) {
		PROTECT(dims = allocVector(VECSXP, 4));
	} else {
		PROTECT(dims = allocVector(VECSXP, 3));
	}
	SET_VECTOR_ELT(dims, 0, mkChar("w"));
	SET_VECTOR_ELT(dims, 1, mkChar("mean"));
	SET_VECTOR_ELT(dims, 2, mkChar("cov"));
	if(hasData) {
		SET_VECTOR_ELT(dims, 3, mkChar("labels"));
	}
	
	vectorToSXP(&w, gslW);
	for(i=0; i<kbis; i++) {
		SET_VECTOR_ELT(cov, i, allocMatrix(REALSXP, d, d));
		SEXP temp = VECTOR_ELT(cov, i);
		matrixToSXP(&temp, gslCov[i]);
	}
	
	SET_VECTOR_ELT(res, 0, w);
	SET_VECTOR_ELT(res, 1, mean);
	SET_VECTOR_ELT(res, 2, cov);
	if(hasData) {
		SET_VECTOR_ELT(res, 3, labs);
	}
	
	setAttrib(res, R_NamesSymbol, dims);	
	
	dropList(selection);
	Free(selection);
	gsl_vector_free(gslW);
	for(i=0; i<kbis; i++) {
		gsl_matrix_free(gslCov[i]);
	}
	Free(gslCov);
	gsl_matrix_free(decomp);
	gsl_permutation_free(perm);
	for(i=0; i<kbis; i++) {
		gsl_matrix_free(inverses[i]);
	}
	Free(inverses);
	if(hasData) {
		UNPROTECT(6);
	} else {
		UNPROTECT(5);
	}
	return(res);
	
}

SEXP buildPlainMod(SEXP list, SEXP inds, SEXP origins) {
	// build a plain model from a list set of models
	// inds will indicate indices of used elements in the list
	// if origins is TRUE : also return a vector indicating to which element in the original list a single component belongs
	
	unsigned i, j;
	// get list length
	int n = length(list);
	int totsize = 0;
	int cursize;
	double *labels;
	
	PROTECT(origins = coerceVector(origins, INTSXP));
	
	// inds is optional : if null, we rebuild inds with 0:n-1
	if(inds == R_NilValue) {
		PROTECT(inds = allocVector(INTSXP, n));
		for(i=0; i<n; i++) {
			INTEGER(inds)[i] = i;
		}
	} else {
		// keep the balance of protects
		PROTECT(inds = coerceVector(inds, INTSXP));
	}
	
	SEXP res;
	if(INTEGER(origins)[0]) {
		PROTECT(res = allocVector(VECSXP, 4));
	} else {
		PROTECT(res = allocVector(VECSXP, 3));
	}
	
	// init to nilvalues
	SET_VECTOR_ELT(res, 0, R_NilValue);
	SET_VECTOR_ELT(res, 1, R_NilValue);
	SET_VECTOR_ELT(res, 2, R_NilValue);
	
	// fill in : we manage pointers to a static structure, no (re)allocation needed.
	// garbage collection ensures pointer vectors that are no more in the tree (one per append operation) are finally dropped.
	for(i=0; i<length(inds); i++) {
		SET_VECTOR_ELT(res, 0, appendRealVector(VECTOR_ELT(res, 0), getListElement(VECTOR_ELT(list, INTEGER(inds)[i]), "w")));
		SET_VECTOR_ELT(res, 1, appendList(VECTOR_ELT(res, 1), getListElement(VECTOR_ELT(list, INTEGER(inds)[i]), "mean")));
		SET_VECTOR_ELT(res, 2, appendList(VECTOR_ELT(res, 2), getListElement(VECTOR_ELT(list, INTEGER(inds)[i]), "cov")));
		cursize = length(getListElement(VECTOR_ELT(list, INTEGER(inds)[i]), "w"));
		if(INTEGER(origins)[0]) {
			if(totsize == 0) {
				labels = Calloc(cursize, double);
				for(j=0; j<cursize; j++) {
					labels[j] = i;
				}
				totsize = cursize;
			} else {
				labels = Realloc(labels, totsize + cursize, double);
				for(j=0; j<cursize; j++) {
					labels[totsize + j] = i;
				}
				totsize += cursize;
			}
		}			
			
	}
	
	// normalize new w vector
	double tot = sum(REAL(VECTOR_ELT(res, 0)), length(VECTOR_ELT(res, 0)), 1);
	double *newval = REAL(VECTOR_ELT(res, 0));
	for(i=0; i<length(VECTOR_ELT(res, 0)); i++) {
		newval[i] = newval[i] / tot;
	}
	
	// create SXP labels if needed
	if(INTEGER(origins)[0]) {
		SET_VECTOR_ELT(res, 3, allocVector(INTSXP, totsize));
		for(i=0; i<totsize; i++) {
			INTEGER(VECTOR_ELT(res, 3))[i] = labels[i];
		}
	}
		
		
	// reset list attributes
	SEXP dims;
	if(INTEGER(origins)[0]) {
		PROTECT(dims = allocVector(VECSXP, 4));
	} else {
		PROTECT(dims = allocVector(VECSXP, 3));
	}
	SET_VECTOR_ELT(dims, 0, mkChar("w"));
	SET_VECTOR_ELT(dims, 1, mkChar("mean"));
	SET_VECTOR_ELT(dims, 2, mkChar("cov"));
	if(INTEGER(origins)[0]) {
		SET_VECTOR_ELT(dims, 3, mkChar("a"));
	}

	setAttrib(res, R_NamesSymbol, dims);
	

	UNPROTECT(4);
	return(res);

}


SEXP appendRealVector(SEXP cont, SEXP cue) {
	// append cue to container	
	SEXP res;
	int l1;
	if(cont == R_NilValue) {
		l1 = 0;
	} else {
		l1 = length(cont);
	}
	int l2 = length(cue);
	int i;
	
	// actual values are used as we assign to *val.
	PROTECT(res = allocVector(REALSXP, l1+l2));
	for(i=0; i<l1; i++) {
		REAL(res)[i] = REAL(cont)[i];
	}
	
	for(i=0; i<l2; i++) {
		REAL(res)[l1+i] = REAL(cue)[i];
	}
	
	UNPROTECT(1);
	return(res);
	
}

SEXP appendList(SEXP cont, SEXP cue) {
	SEXP res;
	// manage initial case
	int l1;
	if(cont == R_NilValue) {
		l1 = 0;
	} else {
		l1 = length(cont);
	}
	int l2 = length(cue);
	int i;
	
	PROTECT(res = allocVector(VECSXP, l1+l2));
	// we do not reallocate vector elements : we assign pointer to a static structure.
	for(i=0; i<l1; i++) {
		SET_VECTOR_ELT(res, i, VECTOR_ELT(cont, i));
	}
	
	for(i=0; i<l2; i++) {
		SET_VECTOR_ELT(res, l1+i, VECTOR_ELT(cue, i));
	}
	
	UNPROTECT(1);
	return(res);	
}




SEXP GSLklut(SEXP mod1, SEXP mod2) {
	// computes klut(mod1 || mod2)
	// for 2 models structured (w, mean, cov)
	// assumes coherent dimensionalities
	int k1 = length(getListElement(mod1, "w"));
	int k2 = length(getListElement(mod2, "w"));
	int d = length(VECTOR_ELT(getListElement(mod1, "mean"), 0));
	int i,j,k;
	
	SEXP inputCov1 = getListElement(mod1, "cov");
	SEXP inputCov2 = getListElement(mod2, "cov");
	//CHANGE : store vector_views to each individual weight, mean, cov
	double *k1coeffs = REAL(getListElement(mod1, "w"));
	double *k2coeffs = REAL(getListElement(mod2, "w"));
	gsl_vector_view k1means[k1];
	gsl_vector_view k2means[k2];
	for(i=0; i<k1; i++) {
		k1means[i] = gsl_vector_view_array(REAL(VECTOR_ELT(getListElement(mod1, "mean"), i)), d);
	}
	for(i=0; i<k2; i++) {
		k2means[i] = gsl_vector_view_array(REAL(VECTOR_ELT(getListElement(mod2, "mean"), i)), d);
	}
	
	// create data structure for the sample
	// NB easier to have a sample per k1 components.
	// rather use Calloc
	
	// CHANGE : store sample by individual line
	gsl_vector ***samp = Calloc(k1, gsl_vector **);
	for(i=0; i<k1; i++) {
		samp[i] = Calloc(2 * d + 1, gsl_vector *);
		for(j=0; j<(2*d+1); j++) {
			samp[i][j] = gsl_vector_alloc(d);
		}
	}
	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(4 * d);
	gsl_matrix *a = gsl_matrix_alloc(d,d);
	gsl_matrix *vecs = gsl_matrix_alloc(d,d);
	gsl_vector *vals = gsl_vector_alloc(d);
	gsl_vector_view view1;
	gsl_matrix_view matview1;
	
	
		
	// perform eigen value decompositions, and fill sample matrix with its results
	for(i=0; i<k1; i++) {
		matview1 = gsl_matrix_view_array(REAL(VECTOR_ELT(inputCov1, i)), d, d);
		gsl_matrix_memcpy(a, &matview1.matrix);
		gsl_eigen_symmv(a, vals, vecs, workspace);
		for(int j=0; j<d; j++) {
			view1 = gsl_matrix_column(vecs, j);
			
			gsl_blas_dscal(sqrt(gsl_vector_get(vals, j)), &view1.vector);
			gsl_blas_dcopy(&k1means[i].vector, samp[i][2*j]);
			gsl_blas_daxpy(1.0, &view1.vector, samp[i][2*j]);			
			
			gsl_blas_dcopy(&k1means[i].vector, samp[i][2*j+1]);
			gsl_blas_daxpy(-1.0, &view1.vector, samp[i][2*j+1]);
		}
		gsl_blas_dcopy(&k1means[i].vector, samp[i][2*d]);
	}
	
	// pre-decomp all covariances
	// data for calculation
	// and precalc dets and inverses
	gsl_matrix **k1inverse = Calloc(k1, gsl_matrix *);
	gsl_matrix **k2inverse = Calloc(k2, gsl_matrix *);
	for(i=0; i<k1; i++) {
		k1inverse[i] = gsl_matrix_alloc(d,d);
	}
	for(i=0; i<k2; i++) {
		k2inverse[i] = gsl_matrix_alloc(d,d);
	}
	double k1dets[k1];
	double k2dets[k2];
	
	gsl_matrix *decomp = gsl_matrix_alloc(d,d);
	gsl_permutation *perm = gsl_permutation_alloc(d);
	int ind;
	
	for(i=0; i<k1; i++) {
		matview1 = gsl_matrix_view_array(REAL(VECTOR_ELT(inputCov1, i)), d, d);
		gsl_matrix_memcpy(decomp, &matview1.matrix);
		gsl_linalg_LU_decomp(decomp,perm, &ind);
		k1dets[i] = gsl_linalg_LU_det(decomp, ind);
		gsl_linalg_LU_invert(decomp, perm, k1inverse[i]);
	}
	
	for(i=0; i<k2; i++) {
		matview1 = gsl_matrix_view_array(REAL(VECTOR_ELT(inputCov2, i)), d, d);
		gsl_matrix_memcpy(decomp, &matview1.matrix);
		gsl_linalg_LU_decomp(decomp,perm, &ind);
		k2dets[i] = gsl_linalg_LU_det(decomp, ind);
		gsl_linalg_LU_invert(decomp, perm, k2inverse[i]);
	}
	
		
	
		
	// calculate KL approximation
	double klapprox = 0.0;
	double constant = pow(2 * M_PI, d/2.0);
	constant = 1/constant;
	double quad;
	double coeff;
	double detcoeff;
	gsl_vector *temp1 = gsl_vector_alloc(d);

	
	for(i=0; i<k1; i++) {
		for(j=0; j<(2*d+1); j++) {
			double curval=0.0;
			for(k=0; k<k1; k++) {
				gsl_blas_daxpy(-1.0, &k1means[k].vector, samp[i][j]);
				gsl_blas_dsymv(CblasUpper, -0.5, k1inverse[k], samp[i][j], 0.0, temp1);
				gsl_blas_ddot(samp[i][j], temp1, &quad);
				quad = exp(quad);
				coeff = k1coeffs[k] * constant;
				detcoeff = 1.0 / pow(k1dets[k], 0.5);
				curval += quad * coeff * detcoeff;
				gsl_blas_daxpy(1.0, &k1means[k].vector, samp[i][j]);
			}
			klapprox += (1.0 / (2.0 * d + 1.0)) * k1coeffs[i] * log(curval);
			curval = 0.0;
			for(k=0; k<k2; k++) {
				gsl_blas_daxpy(-1.0, &k2means[k].vector, samp[i][j]);
				gsl_blas_dsymv(CblasUpper, -0.5, k2inverse[k], samp[i][j], 0.0, temp1);
				gsl_blas_ddot(samp[i][j], temp1, &quad);
				quad = exp(quad);
				coeff = k2coeffs[k] * constant;
				detcoeff = 1.0 / pow(k2dets[k], 0.5);
				curval += quad * coeff * detcoeff;
				gsl_blas_daxpy(1.0, &k2means[k].vector, samp[i][j]);
			}
			klapprox -= (1.0 / (2.0 * d + 1.0)) * k1coeffs[i] * log(curval);
			
		}
	}
			

	// deallocate everything
	for(i=0; i<k1; i++) {
		for(j=0; j<(2*d+1); j++) {
			gsl_vector_free(samp[i][j]);
		}
		Free(samp[i]);
		gsl_matrix_free(k1inverse[i]);
	}
	for(i=0; i<k2; i++) {
		gsl_matrix_free(k2inverse[i]);
	}
	Free(samp);
	Free(k1inverse);
	Free(k2inverse);
	gsl_eigen_symmv_free(workspace);
	gsl_matrix_free(a);
	gsl_matrix_free(vecs);
	gsl_vector_free(vals);
	gsl_vector_free(temp1);
	gsl_matrix_free(decomp);
	gsl_permutation_free(perm);
	
	SEXP res;
	PROTECT(res = allocVector(REALSXP, 1));
	REAL(res)[0] = klapprox;
	UNPROTECT(1);
	return(res);
		
}
		
	
SEXP GSLklutSet(SEXP mod1, SEXP set) {
	// computes KLUT of a model wrt a set of models
	// sample shall be computed once.
	int nbset = length(set);
	int k1 = length(getListElement(mod1, "w"));
	int k2[nbset];
	
	
	int d = length(VECTOR_ELT(getListElement(mod1, "mean"), 0));
	int i,j,k,l;
	
	for(i=0; i<nbset; i++) {
		k2[i] = length(getListElement(VECTOR_ELT(set, i), "w"));
	}
	
	SEXP inputCov1 = getListElement(mod1, "cov");
	SEXP inputCov2[nbset];
	for(i=0; i<nbset; i++) {
		inputCov2[i] = getListElement(VECTOR_ELT(set, i), "cov");
	}
	
	//CHANGE : store vector_views to each individual weight, mean, cov
	double *k1coeffs = REAL(getListElement(mod1, "w"));
	double **k2coeffs;
	k2coeffs = Calloc(nbset, double *);
	for(i=0; i<nbset; i++) {
		k2coeffs[i] = REAL(getListElement(VECTOR_ELT(set, i), "w"));
	}
	
	gsl_vector_view k1means[k1];
	gsl_vector_view **k2means = Calloc(nbset, gsl_vector_view *);
	for(i=0; i<nbset; i++) {
		k2means[i] = Calloc(k2[i], gsl_vector_view);
	}
	
	
	for(i=0; i<k1; i++) {
		k1means[i] = gsl_vector_view_array(REAL(VECTOR_ELT(getListElement(mod1, "mean"), i)), d);
	}
	
	for(i=0; i<nbset; i++) {
		for(j=0; j<k2[i]; j++) {
			k2means[i][j] = gsl_vector_view_array(REAL(VECTOR_ELT(getListElement(VECTOR_ELT(set, i), "mean"), j)), d);
		}
	}
	
	// create data structure for the sample
	// NB easier to have a sample per k1 components.
	// rather use Calloc
	
	// CHANGE : store sample by individual line
	gsl_vector ***samp = Calloc(k1, gsl_vector **);
	for(i=0; i<k1; i++) {
		samp[i] = Calloc(2 * d + 1, gsl_vector *);
		for(j=0; j<(2*d+1); j++) {
			samp[i][j] = gsl_vector_alloc(d);
		}
	}
	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(4 * d);
	gsl_matrix *a = gsl_matrix_alloc(d,d);
	gsl_matrix *vecs = gsl_matrix_alloc(d,d);
	gsl_vector *vals = gsl_vector_alloc(d);
	gsl_vector_view view1;
	gsl_matrix_view matview1;
	
	
		
	// perform eigen value decompositions, and fill sample matrix with its results
	for(i=0; i<k1; i++) {
		matview1 = gsl_matrix_view_array(REAL(VECTOR_ELT(inputCov1, i)), d, d);
		gsl_matrix_memcpy(a, &matview1.matrix);
		gsl_eigen_symmv(a, vals, vecs, workspace);
		for(int j=0; j<d; j++) {
			view1 = gsl_matrix_column(vecs, j);
			
			gsl_blas_dscal(sqrt(gsl_vector_get(vals, j)), &view1.vector);
			gsl_blas_dcopy(&k1means[i].vector, samp[i][2*j]);
			gsl_blas_daxpy(1.0, &view1.vector, samp[i][2*j]);			
			
			gsl_blas_dcopy(&k1means[i].vector, samp[i][2*j+1]);
			gsl_blas_daxpy(-1.0, &view1.vector, samp[i][2*j+1]);
		}
		gsl_blas_dcopy(&k1means[i].vector, samp[i][2*d]);
	}
	
	// pre-decomp all covariances
	// data for calculation
	// and precalc dets and inverses
	gsl_matrix **k1inverse = Calloc(k1, gsl_matrix *);
	
	//gsl_matrix **k2inverse = Calloc(k2, gsl_matrix *);
	gsl_matrix ***k2inverse = Calloc(nbset, gsl_matrix **);
	
	for(i=0; i<k1; i++) {
		k1inverse[i] = gsl_matrix_alloc(d,d);
	}
	for(i=0; i<nbset; i++) {
		k2inverse[i] = Calloc(k2[i], gsl_matrix *);
		for(j=0; j<k2[i]; j++) {
			k2inverse[i][j] = gsl_matrix_alloc(d,d);
		}
	}
	double k1dets[k1];
	//double k2dets[k2];
	double **k2dets;
	k2dets = Calloc(nbset, double *);
	for(i=0; i<nbset; i++) {
		k2dets[i] = Calloc(k2[i], double);
	}
	
	gsl_matrix *decomp = gsl_matrix_alloc(d,d);
	gsl_permutation *perm = gsl_permutation_alloc(d);
	int ind;
	
	for(i=0; i<k1; i++) {
		matview1 = gsl_matrix_view_array(REAL(VECTOR_ELT(inputCov1, i)), d, d);
		gsl_matrix_memcpy(decomp, &matview1.matrix);
		gsl_linalg_LU_decomp(decomp,perm, &ind);
		k1dets[i] = gsl_linalg_LU_det(decomp, ind);
		gsl_linalg_LU_invert(decomp, perm, k1inverse[i]);
	}
	
	for(i=0; i<nbset; i++) {
		for(j=0; j<k2[i]; j++) {
			matview1 = gsl_matrix_view_array(REAL(VECTOR_ELT(inputCov2[i], j)), d, d);
			gsl_matrix_memcpy(decomp, &matview1.matrix);
			gsl_linalg_LU_decomp(decomp,perm, &ind);
			k2dets[i][j] = gsl_linalg_LU_det(decomp, ind);
			gsl_linalg_LU_invert(decomp, perm, k2inverse[i][j]);
		}
	}
	
		
	
		
	// calculate KL approximation
	//double klapprox = 0.0;
	double klapprox[nbset];
	for(i=0; i<nbset; i++) klapprox[i] = 0.0;
	
	double constant = pow(2 * M_PI, d/2.0);
	constant = 1/constant;
	double quad;
	double coeff;
	double detcoeff;
	gsl_vector *temp1 = gsl_vector_alloc(d);

	
	for(i=0; i<k1; i++) {
		for(j=0; j<(2*d+1); j++) {
			double curval=0.0;
			for(k=0; k<k1; k++) {
				gsl_blas_daxpy(-1.0, &k1means[k].vector, samp[i][j]);
				gsl_blas_dsymv(CblasUpper, -0.5, k1inverse[k], samp[i][j], 0.0, temp1);
				gsl_blas_ddot(samp[i][j], temp1, &quad);
				quad = exp(quad);
				coeff = k1coeffs[k] * constant;
				detcoeff = 1.0 / pow(k1dets[k], 0.5);
				curval += quad * coeff * detcoeff;
				gsl_blas_daxpy(1.0, &k1means[k].vector, samp[i][j]);
			}
			for(k=0; k<nbset; k++) {
				klapprox[k] += (1.0 / (2.0 * d + 1.0)) * k1coeffs[i] * log(curval);
			}
			for(l=0; l<nbset; l++) {
				curval = 0.0;
				for(k=0; k<k2[l]; k++) {
					gsl_blas_daxpy(-1.0, &k2means[l][k].vector, samp[i][j]);
					gsl_blas_dsymv(CblasUpper, -0.5, k2inverse[l][k], samp[i][j], 0.0, temp1);
					gsl_blas_ddot(samp[i][j], temp1, &quad);
					quad = exp(quad);
					coeff = k2coeffs[l][k] * constant;
					detcoeff = 1.0 / pow(k2dets[l][k], 0.5);
					curval += quad * coeff * detcoeff;
					gsl_blas_daxpy(1.0, &k2means[l][k].vector, samp[i][j]);
				}
				klapprox[l] -= (1.0 / (2.0 * d + 1.0)) * k1coeffs[i] * log(curval);
			}
			
		}
	}
			

	// deallocate everything
	Free(k2coeffs);
	for(i=0; i<nbset; i++) {
		Free(k2means[i]);
	}
	Free(k2means);
	

	
	for(i=0; i<nbset; i++) {
		Free(k2dets[i]);
	}
	Free(k2dets);
	
	for(i=0; i<k1; i++) {
		for(j=0; j<(2*d+1); j++) {
			gsl_vector_free(samp[i][j]);
		}
		Free(samp[i]);
		gsl_matrix_free(k1inverse[i]);
	}
	for(i=0; i<nbset; i++) {
		for(j=0; j<k2[i]; j++) {
			gsl_matrix_free(k2inverse[i][j]);
		}
		Free(k2inverse[i]);
	}
	Free(k2inverse);
	Free(samp);
	Free(k1inverse);
	gsl_eigen_symmv_free(workspace);
	gsl_matrix_free(a);
	gsl_matrix_free(vecs);
	gsl_vector_free(vals);
	gsl_vector_free(temp1);
	gsl_matrix_free(decomp);
	gsl_permutation_free(perm);
	
	SEXP res;
	PROTECT(res = allocVector(REALSXP, nbset));
	for(i=0; i<nbset; i++) {
		REAL(res)[i] = klapprox[i];
	}
	UNPROTECT(1);
	return(res);


}	
	

void printGSLvector(gsl_vector *vec) {
	int sz = vec->size;
	int i;
	for(i=0; i<sz; i++) {
		Rprintf("%e ", gsl_vector_get(vec, i));
	}
	Rprintf("\n");
}

void printGSLmatrix(gsl_matrix *mat) {
	int sz1 = mat->size1;
	int sz2 = mat->size2;
	int i,j;
	for(i=0; i<sz1; i++) {
		for(j=0; j<sz2; j++) {
			Rprintf("%e ", gsl_matrix_get(mat, i, j));
		}
		Rprintf("\n");
	}
}
		
		
void printDoubleVector(double *tab, int size) {
	int i;
	for(i=0; i<size; i++) {
		Rprintf("%f ", tab[i]);
	}
	Rprintf("\n");
}	
		
SEXP getTimestamp() {
	// get current timestamp in format (sec, microsec)
	struct timeval *tv;
	tv = (struct timeval *)malloc(sizeof(struct timeval));
	
	gettimeofday(tv, NULL);
	
	SEXP R_res;
	PROTECT(R_res = allocVector(INTSXP, 2));
	INTEGER(R_res)[0] = tv->tv_sec;
	INTEGER(R_res)[1] = tv->tv_usec;
	
	UNPROTECT(1);
	free(tv);
	return(R_res);
}


SEXP getElapsed(SEXP stamp) {
	// get time elapsed since stamp in parameter
	// return result in double valued seconds
	
	struct timeval *tv;
	tv = (struct timeval *)malloc(sizeof(struct timeval));
	
	gettimeofday(tv, NULL);
	
	tv->tv_sec -= INTEGER(stamp)[0];
	if(tv->tv_usec >= INTEGER(stamp)[1]) {
		tv->tv_usec -= INTEGER(stamp)[1];
	} else {
		tv->tv_usec = 1000000 + tv->tv_usec - INTEGER(stamp)[1];
		tv->tv_sec--;
	}	
	
	double res = tv->tv_sec + tv->tv_usec / 1000000.0;
	SEXP R_res;
	PROTECT(R_res = allocVector(REALSXP, 1));
	REAL(R_res)[0] = res;
	
	UNPROTECT(1);
	free(tv);
	
	
	
	return(R_res);
}


double sum(double *tab, int size, int stride) {
	int i;
	double tot=0.0;
	for(i=0; i<size; i++) {
		tot += tab[i*stride];
	}
	return(tot);
}


SEXP pointwise(SEXP mat1, SEXP mat2) {
	// multiplication of the vectors in mat1 and mat2
	// mat1 has r1 rows and c1 cols, respectively for mat2
	// must have : r1=c2 and r2=c1
	// we multiply the r1 rows of mat1 pointwise against the c2 cols of mat2
	// thus r1 values are obtained as output.
	
	int r1 = INTEGER(getAttrib(mat1, R_DimSymbol))[0];
	int c1 = INTEGER(getAttrib(mat1, R_DimSymbol))[1];
	
	PROTECT(mat1=coerceVector(mat1, REALSXP));
	PROTECT(mat2=coerceVector(mat2, REALSXP));
	
	int i;
	double vals[r1];
	
	gsl_matrix_view matview1;
	gsl_matrix_view matview2;
	gsl_vector_view view1;
	gsl_vector_view view2;
	
	// we read the matrixes "reverse" as storage is inversed in R
	matview1 = gsl_matrix_view_array(REAL(mat1), c1, r1);
	matview2 = gsl_matrix_view_array(REAL(mat2), r1, c1);
	
	for(i=0; i<r1; i++) {
		view1 = gsl_matrix_column(&matview1.matrix, i);
		view2 = gsl_matrix_row(&matview2.matrix, i);
		gsl_blas_ddot(&view1.vector, &view2.vector, &(vals[i]));
	}
	
	
	SEXP res;
	PROTECT(res=allocVector(REALSXP, r1));
	for(i=0; i<r1; i++) {
		REAL(res)[i] = vals[i];
	}
	
	UNPROTECT(3);
	return(res);
}	
	
	
int gramschmidt(gsl_matrix *mat) {
	// perform gramschmidt diagonalisation
	int i,j;
	int d = mat->size2;
	int n = mat->size1;
	double norm1;
	
	gsl_vector_view view1;
	gsl_vector_view view2;
	
	for(i=0; i<d; i++) {
		view1 = gsl_matrix_column(mat, i);
		for(j=0; j<i; j++) {
			view2 = gsl_matrix_column(mat, j);
			gsl_blas_ddot(&view1.vector, &view2.vector, &norm1);
			gsl_blas_daxpy(-norm1, &view2.vector, &view1.vector);
		}
		gsl_blas_ddot(&view1.vector, &view1.vector, &norm1);
		norm1 = sqrt(norm1);
		gsl_blas_dscal(1.0/norm1, &view1.vector);
	}
	
	return(0);
}



SEXP R_gramschmidt(const SEXP R_mat) {
	// use SXPtoMatrix
	// gramschmidt modifies original values : no pb, because gslmatrix locally created.
	int nrow = INTEGER(getAttrib(R_mat, R_DimSymbol))[0];
	int ncol = INTEGER(getAttrib(R_mat, R_DimSymbol))[1];	
	
	SEXP R_out;
	PROTECT(R_out = allocMatrix(REALSXP, nrow, ncol));
	
	gsl_matrix *mat = gsl_matrix_alloc(nrow, ncol);
	SXPtoMatrix(mat, R_mat);
	
	gramschmidt(mat);
	
	matrixToSXP(&R_out, mat);
	
	UNPROTECT(1);
	gsl_matrix_free(mat);
	return(R_out);
}


void upperComplete(gsl_matrix *mat) {
	// complete lower triangular part in a matrix(assumed to be symmetric)
	unsigned i,j;
	double cur;
	int dim = mat->size1;
	for(i=0; i<(dim-1); i++) {
		for(j=(i+1); j<dim; j++) {

			cur = gsl_matrix_get(mat, i, j);
			gsl_matrix_set(mat, j, i, cur);
		}
	}
}

void lowerComplete(double *mat, int n) {
	// complete upper triangle
	for(int i=0; i<n; i++) {
		for(int j=0; j<i; j++) {
			mat[i*n+j] = mat[j*n+i];
		}
	}
}

void fillUpper(double *mat, int d) {
	for(int i=1; i<d; i++) {
		for(int j=0; j<i; j++) {
			mat[j+d*i] = mat[i+d*j];
		}
	}
}

void fillLower(double *mat, int d) {
	for(int i=1; i<d; i++) {
		for(int j=0; j<i; j++) {
			mat[i+d*j] = mat[j+d*i];
		}
	}
}

SEXP sort_index(SEXP vec, SEXP order) {
	double *ptr = REAL(coerceVector(vec, REALSXP));
	int len = length(vec);
	SEXP R_res;
	int *ptr2;
	int i, c_order;

	// 0 for asc, 1 for desc
	c_order = INTEGER(coerceVector(order, INTSXP))[0];
	
	gsl_vector_view gslvec = gsl_vector_view_array(ptr, len);
	gsl_permutation *res = gsl_permutation_alloc(len);
	
	gsl_sort_vector_index(res, &gslvec.vector);
	
	PROTECT(R_res = allocVector(INTSXP, len));
	ptr2 = INTEGER(R_res);

	for(i=0; i<len; i++) {
		// we return indices from 1 to n (instead of 0 to n-1)
		if(c_order) {
			ptr2[len-1-i] = res->data[i] + 1;
		} else {
			ptr2[i] = res->data[i] + 1;
		}
	}
	
	gsl_permutation_free(res);
	UNPROTECT(1);
	return(R_res);
}

// sample sizeof(samp) elements from 1..n
void GSLsample(int n, int size, int *samp) {
	int i, res, ok;
	double cur, lim;
	GetRNGstate();
	for(i=0; i<size; i++) {
		ok=0;
		while(!ok) {
			cur = runif(0.0,1.0);
			lim = 1.0/(double)n;
			res = 1;

			while(cur > lim) {
				lim += 1.0/(double)n;
				res++;
			}
			if(!contains(res, size, samp)) ok = 1;
		}
		samp[i] = res;
	}
	PutRNGstate();
}

// sample samp elements from 1..n, return list.
SEXP sample(SEXP n, SEXP size) {
	int c_n = INTEGER(coerceVector(n, INTSXP))[0];
	int c_size = INTEGER(coerceVector(size, INTSXP))[0];
	int i;
	int *samp = Calloc(c_size, int);
	GSLsample(c_n, c_size, samp);
	
	SEXP res;
	PROTECT(res=allocVector(INTSXP, c_size));
	for(i=0; i<c_size; i++) {
		INTEGER(res)[i] = samp[i];
	}
	
	Free(samp);
	UNPROTECT(1);
	return(res);
}

int contains(int seekval, int size, int *vals) {
	int i;
	for(i=0; i<size; i++) {
		if(vals[i] == seekval) return 1;
	}
	return 0;
}



SEXP control(SEXP points, SEXP dmin) {
	int d = INTEGER(getAttrib(points, R_DimSymbol))[1];
	int n = INTEGER(getAttrib(points, R_DimSymbol))[0];
	int i,j;
	double Cdmin = INTEGER(coerceVector(dmin, INTSXP))[0];

	gsl_matrix *GSLpoints = gsl_matrix_alloc(n, d);
	for(i=0; i<n; i++) {
		for(j=0; j<d; j++) {
			gsl_matrix_set(GSLpoints, i, j, REAL(points)[j*n + i]);
		}
	}
	
	gsl_vector_view curpoint1;
	gsl_vector_view curpoint2;
	
	SEXP indic;
	PROTECT(indic=allocVector(INTSXP, n));
	
	int isFarFromOthers;
	

	for(i=0; i<n; i++) {
		curpoint1 = gsl_matrix_row(GSLpoints, i);
		isFarFromOthers = 1;
		for(j=0; j<n; j++) {
			curpoint2 = gsl_matrix_row(GSLpoints, j);
			if(i != j) {
				if(getDistance(&curpoint1.vector, &curpoint2.vector) < Cdmin) isFarFromOthers = 0;
			}
		}
		INTEGER(indic)[i] = isFarFromOthers;
	}
	
	gsl_matrix_free(GSLpoints);
	UNPROTECT(1);
	return(indic);
	
}



double getDistance(gsl_vector *vec1, gsl_vector *vec2) {
	double total=0.0;
	
	for(int i=0; i<vec1->size; i++) {
		total += pow(gsl_vector_get(vec1, i) - gsl_vector_get(vec2, i), 2.0);
	}
	
	return(sqrt(total));
	
}


void getColumnNorms(gsl_matrix *mat, gsl_vector *vec) {
	int d = mat->size2;
	int n = mat->size1;
	int i,j;
	double cur;

	for(i=0; i<d; i++) {
		cur=0.0;
		for(j=0; j<n; j++) {
			cur += pow(gsl_matrix_get(mat, j, i), 2.0);
		}
		gsl_vector_set(vec, i, sqrt(cur));
	}
}

void getCovariance(gsl_matrix *elems, gsl_matrix *cov) {
	// get sample covariance of elements in lines of elems
	// store in cov
	int i;
	gsl_vector *temp_mean = gsl_vector_alloc(elems->size2);
	gsl_vector *temp_acc = gsl_vector_alloc(elems->size2);
	getMean(elems, temp_mean);
	gsl_vector_view view1;
	
	gsl_matrix_set_zero(cov);
	
	for(i=0; i<elems->size1; i++) {
		view1 = gsl_matrix_row(elems, i);
		gsl_vector_memcpy(temp_acc, &view1.vector);
		gsl_blas_daxpy(-1.0, temp_mean, temp_acc);
		gsl_blas_dsyr(CblasUpper, 1.0 / (double)elems->size1, temp_acc, cov);
	}
	
}

void getMean(gsl_matrix *elems, gsl_vector *mean) {
	// get mean of lines in elems
	int i;
	gsl_vector_set_zero(mean);
	gsl_vector_view view1;
	
	for(i=0; i<elems->size1; i++) {
		view1 = gsl_matrix_row(elems, i);
		gsl_blas_daxpy(1.0, &view1.vector, mean);
	}
	gsl_vector_scale(mean, 1.0 / (double)elems->size1);
}


SEXP rDirichlet(SEXP K, SEXP R_alpha) {
	int k = INTEGER(coerceVector(K, INTSXP))[0];
	
	gsl_rng_env_setup();
	gsl_rng_default_seed = INTEGER(getTimestamp())[1];
	
	const gsl_rng_type * T;
	T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(T);
	
	double alpha[k];
	for(int i=0; i<k; i++) alpha[i] = REAL(coerceVector(R_alpha, REALSXP))[0];
	double theta[k];
	
	gsl_ran_dirichlet(r, k, alpha, theta);
	
	SEXP res;
	PROTECT(res = allocVector(REALSXP, k));
	for(int i=0; i<k; i++) REAL(res)[i] = theta[i];
	
	UNPROTECT(1);
	gsl_rng_free(r);
	return(res);
}

// dist between 2 groups of elements, wrt to a mahalanobis metric
SEXP gdist(SEXP g1, SEXP g2, SEXP metric) {
	PROTECT(g1 = coerceVector(g1, REALSXP));
	PROTECT(g2 = coerceVector(g2, REALSXP));
	// metric is a list:
	// if length=1, dists between elts in g1 and g2 are computed wrt this metric,
	// if length=length(g2), metrics are associated with elements in g2.
	//PROTECT(metric = coerceVector(metric, REALSXP));
	SEXP curmetric;
	double *c_metric;
	int metricMode = length(metric);

	int n1 = INTEGER(getAttrib(g1, R_DimSymbol))[0];
	int n2 = INTEGER(getAttrib(g2, R_DimSymbol))[0];
	int d = INTEGER(getAttrib(g1, R_DimSymbol))[1];	
	double *c_g1 = REAL(g1);
	double *c_g2 = REAL(g2);

	//double *c_metric = REAL(metric);
	double *vec1 = calloc(d, sizeof(double));
	double *vec2 = calloc(d, sizeof(double));
	double *c_inverse = calloc(d*d, sizeof(double));
	int i1=1;
	double alpha;
	double beta=0.0;
	double lndet = 0.0;

	SEXP dists;
	PROTECT(dists=allocMatrix(REALSXP,n1,n2));
	double *c_dists=REAL(dists);

	if(metricMode == 1) {
		c_metric = REAL(VECTOR_ELT(metric, 0));
	}

	for(int j=0; j<n2; j++) {
		// inverted loop so that matrices are decomposed only once
		if(metricMode == n2) {
			c_metric = REAL(VECTOR_ELT(metric, j));
			symdecomp(c_metric, c_inverse, &lndet, d, "general");
		}
		for(int i=0; i<n1; i++) {
			F77_CALL(dcopy)(&d, c_g1+i, &n1, vec1, &i1);
			alpha=-1.0;
			F77_CALL(daxpy)(&d, &alpha, c_g2+j, &n2, vec1, &i1);
			alpha=1.0;

			F77_CALL(dsymv)("L", &d, &alpha, c_inverse, &d, vec1, &i1, &beta, vec2, &i1);
			c_dists[j*n1+i] = sqrt(F77_CALL(ddot)(&d, vec1, &i1, vec2, &i1));
		}
	}

	free(vec1);
	free(vec2);
	UNPROTECT(3);
	return(dists);

}


