// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include "mppca.h"

void writeMatrix(gsl_matrix *mat) {
	QFile file("/Users/Pierrick/R/log.txt");
	file.open(QIODevice::Append);
	
	QTextStream out(&file);	
	// write line to identify output
	out << "logged at" << QDateTime::currentDateTime().toString() << "\n";
	out << "-------------------\n";
	
	int i,j;
	QString curnum;
	for(i=0; i<mat->size1; i++) {
		for(j=0; j<mat->size2; j++) {
			curnum.setNum(gsl_matrix_get(mat, i, j), 'f');
			out << curnum << " ";
		}
		out << "\n";
	}
	out << "\n";
}


// HOWTO get rows from a R matrix into a GSL program :
// - invert usage of row and col count to define data matrix
// - get vectors from gsl_matrix_column instead of row.
// Logically, data is transposed in our program.

void alloc(SEXP data, SEXP ncomp, SEXP qmax) {
	K = INTEGER(coerceVector(ncomp, INTSXP))[0];
	n = INTEGER(getAttrib(data, R_DimSymbol))[0];
	d = INTEGER(getAttrib(data, R_DimSymbol))[1];

	double *dataptr = REAL(coerceVector(data, REALSXP));
	dataview = gsl_matrix_view_array(dataptr, d, n);
	
	int c_qmax = INTEGER(coerceVector(qmax, INTSXP))[0];
	q = ((d-1)<c_qmax) ? d-1 : c_qmax;
	
	stat_z = gsl_matrix_alloc(n,K);
	stat_norms = gsl_vector_alloc(n);
	stat_lb = gsl_vector_alloc(8);
	stat_trace = gsl_matrix_alloc(d,q);
	stat_nk = gsl_vector_calloc(K);
	stat_yk = gsl_matrix_calloc(K, d);
	stat_sk = gsl_matrix_calloc(K, q);
	stat_syk = Calloc(K, gsl_matrix *);
	stat_Sk = Calloc(K, gsl_matrix *);
	for(i=0; i<K; i++) {
		stat_syk[i] = gsl_matrix_calloc(d,q);
		stat_Sk[i] = gsl_matrix_calloc(q,q);
	}
	stat_yk2 = gsl_vector_calloc(K);
	
	param_mustar = gsl_matrix_alloc(K, d);
	model_alpha = gsl_vector_alloc(K);
	model_nua = gsl_vector_alloc(K);
	model_nub = gsl_matrix_alloc(K,q);
	model_numoment = gsl_matrix_alloc(K, q);
	model_taua = gsl_vector_alloc(K);
	model_taub = gsl_vector_alloc(K);
	model_taumoment = gsl_vector_alloc(K);
	model_mumean = gsl_matrix_alloc(K,d);
	model_musigma = Calloc(K, gsl_vector *);
	model_wmean = Calloc(K, gsl_matrix *);
	model_xsigma = Calloc(K, gsl_matrix *);
	model_xmean = Calloc(K, gsl_matrix *);
	model_wsigma = Calloc(K, gsl_matrix *);	
	
	for(i=0; i<K; i++) {
		model_musigma[i] = gsl_vector_alloc(d);
		model_wsigma[i] = gsl_matrix_alloc(q,q);
		model_wmean[i] = gsl_matrix_alloc(d,q);
		model_xsigma[i] = gsl_matrix_alloc(q,q);
		model_xmean[i] = gsl_matrix_alloc(n,q);
	}
	
	temp_mat1 = gsl_matrix_alloc(d,d);
	temp_mat2 = gsl_matrix_alloc(d,d);
	temp_mat3 = gsl_matrix_alloc(q,q);
	temp_mat4 = gsl_matrix_alloc(q,q);
	temp_mat5 = gsl_matrix_alloc(q,q);
	temp_mat6 = gsl_matrix_alloc(d,q);
	temp_mat7 = gsl_matrix_alloc(q,q);
	temp_perm1 = gsl_permutation_alloc(d);
	temp_perm3 = gsl_permutation_alloc(q);
	temp_perm4 = gsl_permutation_alloc(q);
	temp_perm5 = gsl_permutation_alloc(n);
	temp_vec1 = gsl_vector_alloc(d);
	temp_vec2 = gsl_vector_alloc(q);
	temp_vec3 = gsl_vector_alloc(q);
	temp_vec4 = gsl_vector_alloc(d);
	temp_vec5 = gsl_vector_alloc(n);
	temp_symmv1 = gsl_eigen_symmv_alloc(q);
	
	mins = Calloc(d, double);
	maxs = Calloc(d, double);
}



void init() {

	for(i=0; i<d; i++) {
		varview1 = gsl_matrix_row(&dataview.matrix, i);
		mins[i] = gsl_vector_min(&varview1.vector);
		maxs[i] = gsl_vector_max(&varview1.vector);
	}

	gsl_matrix_set_zero(stat_z);
	param_alpha = 0.001;
	param_nua = 0.001;
	param_nub = 0.001;
	param_numoment = 1.0;
	param_taua = 0.001;
	param_taub = 0.001;
	param_taumoment = 1.0;
	param_nustar = 0.001;	
	
	GetRNGstate();
	for(i=0; i<K; i++) {
		for(j=0; j<d; j++) {
			val1 = 0.0;
			gsl_matrix_set(param_mustar, i, j, runif(mins[j] - val1, maxs[j] + val1));
		}
	}
	PutRNGstate();	
	
	gsl_vector_set_all(model_alpha, param_alpha);
	
	
	
	gsl_vector_set_all(model_nua, param_nua);
	gsl_matrix_set_all(model_nub, param_nub);
	gsl_matrix_set_all(model_numoment, param_numoment);
	
	gsl_vector_set_all(model_taua, param_taua);
	gsl_vector_set_all(model_taub, param_taub);
	gsl_vector_set_all(model_taumoment, param_taumoment);
	
	
	gsl_matrix_memcpy(model_mumean, param_mustar);
	for(i=0; i<K; i++) {
		gsl_vector_set_zero(model_musigma[i]);

		gsl_vector_set_all(model_musigma[i], 1.0/param_nustar);
	}

	GetRNGstate();
	
	// init all comps to the same basis
	for(i=0; i<q; i++) {
		for(j=0; j<d; j++) {
			gsl_matrix_set(temp_mat6, j, i, rnorm(0.0, 1.0));
		}
	}

	gramschmidt(temp_mat6);
	
		
	for(i=0; i<K; i++) {
		
		gsl_matrix_set_zero(model_wsigma[i]);
		facview1 = gsl_matrix_diagonal(model_wsigma[i]);
		gsl_vector_set_all(&facview1.vector, 1.0);
		
		gsl_matrix_set_zero(model_xsigma[i]);
		facview1 = gsl_matrix_diagonal(model_xsigma[i]);

		gsl_vector_set_all(&facview1.vector, 0.0);

				
		
		for(j=0; j<q; j++) {
			for(k=0; k<d; k++) {
				gsl_matrix_set(model_wmean[i], k, j, rnorm(0.0, 1.0));
			}
		}
		
		
		gramschmidt(model_wmean[i]);
		
		
		for(j=0; j<q; j++) {			
			for(k=0; k<n; k++) {
				gsl_matrix_set(model_xmean[i], k, j, rnorm(0.0, 1.0));
			}
		}

		

		gsl_matrix_set_all(model_xmean[i], 0.0);
		
	}
	PutRNGstate();
	
	
	// init norms
	for(i=0; i<n; i++) {
		indview1 = gsl_matrix_column(&dataview.matrix, i);
		gsl_blas_ddot(&indview1.vector, &indview1.vector, &val1);
		gsl_vector_set(stat_norms, i, val1);
	}
	
}


void initWithPrior(SEXP xmean, SEXP wmean, SEXP mumean) {
	// replace model_xmean and model_wmean
	// dimensions are assumed to be consistent w.r.t input data
	double *xptr; 
	double *wptr;
	double *muptr;
	
	for(i=0; i<K; i++) {
		xptr = REAL(coerceVector(VECTOR_ELT(xmean, i), REALSXP));
		wptr = REAL(coerceVector(VECTOR_ELT(wmean, i), REALSXP));
		muptr = REAL(coerceVector(VECTOR_ELT(mumean, i), REALSXP));
		
		gsl_matrix_view xview = gsl_matrix_view_array(xptr, q, n);
		gsl_matrix_view wview = gsl_matrix_view_array(wptr, q, d);
		gsl_vector_view muview = gsl_vector_view_array(muptr, d);

		indview1 = gsl_matrix_row(model_mumean, i);
		gsl_vector_memcpy(&indview1.vector, &muview.vector);
		indview1 = gsl_matrix_row(param_mustar, i);
		gsl_vector_memcpy(&indview1.vector, &muview.vector);
				

		for(j=0; j<q; j++) {
			varview1 = gsl_matrix_row(&xview.matrix, j);
			varview2 = gsl_matrix_column(model_xmean[i], j);
			gsl_vector_memcpy(&varview2.vector, &varview1.vector);

			indview1 = gsl_matrix_row(&wview.matrix, j);
			indview2 = gsl_matrix_column(model_wmean[i], j);
			gsl_vector_memcpy(&indview2.vector, &indview1.vector);
		}
		
	}

}

/*
void initWithElements() {
	int elements[K];
	GSLsample(n, K, elements);
	
	for(i=0; i<K; i++) {
		indview1 = gsl_matrix_row(param_mustar, i);
		indview2 = gsl_matrix_column(&dataview.matrix, elements[i] - 1);
		gsl_vector_memcpy(&indview1.vector, &indview2.vector);
	}
	gsl_matrix_memcpy(model_mumean, param_mustar);
}
*/

void initWithReadyModel(SEXP inmod) {
	double *ptr;
	
	
	for(i=0; i<K; i++) {
		ptr = REAL(coerceVector(getListElement(inmod, "alpha"), REALSXP));
		gsl_vector_set(model_alpha, i, ptr[i]);

		
		ptr = REAL(coerceVector(getListElement(inmod, "nua"), REALSXP));
		gsl_vector_set(model_nua, i, ptr[i]);
		
		ptr = REAL(coerceVector(getListElement(inmod, "taumoment"), REALSXP));
		gsl_vector_set(model_taumoment, i, ptr[i]);
		
		ptr = REAL(coerceVector(getListElement(inmod, "taua"), REALSXP));
		gsl_vector_set(model_taua, i, ptr[i]);
		
		ptr = REAL(coerceVector(getListElement(inmod, "taub"), REALSXP));
		gsl_vector_set(model_taub, i, ptr[i]);
		
		ptr = REAL(coerceVector(VECTOR_ELT(getListElement(inmod, "nub"), i), REALSXP));
		facview1 = gsl_vector_view_array(ptr, q);
		facview2 = gsl_matrix_row(model_nub, i);
		gsl_vector_memcpy(&facview2.vector, &facview1.vector);
		
		ptr = REAL(coerceVector(VECTOR_ELT(getListElement(inmod, "numoment"), i), REALSXP));
		facview1 = gsl_vector_view_array(ptr, q);
		facview2 = gsl_matrix_row(model_numoment, i);
		gsl_vector_memcpy(&facview2.vector, &facview1.vector);
		
		ptr = REAL(coerceVector(VECTOR_ELT(getListElement(inmod, "mumean"), i), REALSXP));
		facview1 = gsl_vector_view_array(ptr, d);

		facview2 = gsl_matrix_row(model_mumean, i);
		gsl_vector_memcpy(&facview2.vector, &facview1.vector);		
		
		ptr = REAL(coerceVector(VECTOR_ELT(getListElement(inmod, "musigma"), i), REALSXP));

		indview1 = gsl_vector_view_array(ptr, d);
		gsl_vector_memcpy(model_musigma[i], &indview1.vector);
		
		
		ptr = REAL(coerceVector(VECTOR_ELT(getListElement(inmod, "mustar"), i), REALSXP));
		indview1 = gsl_vector_view_array(ptr, d);
		indview2 = gsl_matrix_row(param_mustar, i);
		gsl_vector_memcpy(&indview2.vector, &indview1.vector);
		
		
		
		ptr = REAL(coerceVector(VECTOR_ELT(getListElement(inmod, "wmean"), i), REALSXP));
		tempview1 = gsl_matrix_view_array(ptr, q, d);
		for(j=0; j<d; j++) {
			facview1 = gsl_matrix_column(&tempview1.matrix, j);
			facview2 = gsl_matrix_row(model_wmean[i], j);
			gsl_vector_memcpy(&facview2.vector, &facview1.vector);
		}

		ptr = REAL(coerceVector(VECTOR_ELT(getListElement(inmod, "wsigma"), i), REALSXP));
		tempview1 = gsl_matrix_view_array(ptr, q, q);
		for(j=0; j<q; j++) {
			facview1 = gsl_matrix_column(&tempview1.matrix, j);
			facview2 = gsl_matrix_row(model_wsigma[i], j);
			gsl_vector_memcpy(&facview2.vector, &facview1.vector);
		}

		
	}
	updateX();
}


void initWithReadyW(SEXP inw) {
	double *ptr;
	ptr = REAL(coerceVector(inw, REALSXP));
	tempview1 = gsl_matrix_view_array(ptr, q, d);
	for(i=0; i<K; i++) {
		for(j=0; j<d; j++) {
			facview1 = gsl_matrix_column(&tempview1.matrix, j);
			facview2 = gsl_matrix_row(model_wmean[i], j);
			gsl_vector_memcpy(&facview2.vector, &facview1.vector);
		}
	}
	// update wsigma to reflect the distribution
	gsl_vector_set_zero(temp_vec2);
	gsl_matrix_set_zero(temp_mat3);
	for(i=0; i<d; i++) {
		facview1 = gsl_matrix_row(model_wmean[0], i);
		gsl_vector_add(temp_vec2, &facview1.vector);
	}
	gsl_vector_scale(temp_vec2, 1.0/(double)d);
	
	for(i=0; i<d; i++) {
		facview1 = gsl_matrix_row(model_wmean[0], i);
		gsl_vector_memcpy(temp_vec3, &facview1.vector);
		gsl_vector_sub(temp_vec3, temp_vec2);
		gsl_blas_dsyr(CblasUpper, 1.0, temp_vec3, temp_mat3);
	}
	upperComplete(temp_mat3);
	gsl_matrix_scale(temp_mat3, 1.0/(double)d);
	for(i=0; i<K; i++) {
		gsl_matrix_memcpy(model_wsigma[i], temp_mat3);
	}
}


void end() {
	// deallocate heap variables and model
	gsl_matrix_free(stat_z);
	gsl_vector_free(stat_norms);
	gsl_vector_free(stat_lb);
	gsl_matrix_free(stat_trace);
	gsl_vector_free(stat_nk);
	gsl_matrix_free(stat_yk);
	gsl_matrix_free(stat_sk);
	for(i=0; i<K; i++) {
		gsl_matrix_free(stat_syk[i]);
		gsl_matrix_free(stat_Sk[i]);
	}
	Free(stat_syk);
	Free(stat_Sk);
	gsl_vector_free(stat_yk2);
	
	gsl_matrix_free(param_mustar);
	gsl_vector_free(model_alpha);
	gsl_vector_free(model_nua);
	gsl_matrix_free(model_nub);
	gsl_matrix_free(model_numoment);
	gsl_vector_free(model_taua);
	gsl_vector_free(model_taub);
	gsl_vector_free(model_taumoment);
	gsl_matrix_free(model_mumean);
	
	for(i=0; i<K; i++) {
		gsl_vector_free(model_musigma[i]);
		gsl_matrix_free(model_wsigma[i]);
		gsl_matrix_free(model_wmean[i]);
		gsl_matrix_free(model_xsigma[i]);
		gsl_matrix_free(model_xmean[i]);
	}
	
	Free(model_musigma);
	Free(model_wsigma);
	Free(model_wmean);
	Free(model_xsigma);
	Free(model_xmean);
	
	gsl_matrix_free(temp_mat1);
	gsl_matrix_free(temp_mat2);
	gsl_matrix_free(temp_mat3);
	gsl_matrix_free(temp_mat4);
	gsl_matrix_free(temp_mat5);
	gsl_matrix_free(temp_mat6);
	gsl_matrix_free(temp_mat7);
	gsl_vector_free(temp_vec1);
	gsl_vector_free(temp_vec2);
	gsl_vector_free(temp_vec3);
	gsl_vector_free(temp_vec4);
	gsl_permutation_free(temp_perm1);
	gsl_permutation_free(temp_perm3);
	gsl_permutation_free(temp_perm4);
	gsl_eigen_symmv_free(temp_symmv1);
	
	Free(mins);
	Free(maxs);
	
}



void updateX() {
	// UPDATE X
	for(i=0; i<K; i++) {
		// COVARIANCE UPDATE
		// reset to unit diag
		gsl_matrix_set_zero(model_xsigma[i]);
		facview1 = gsl_matrix_diagonal(model_xsigma[i]);
		gsl_vector_set_all(&facview1.vector, 1.0);
		
		gsl_matrix_memcpy(temp_mat3, model_wsigma[i]);
		gsl_matrix_scale(temp_mat3, gsl_vector_get(model_taumoment, i) * d);
		gsl_matrix_add(model_xsigma[i], temp_mat3);
		for(j=0; j<d; j++) {

			facview1 = gsl_matrix_row(model_wmean[i], j);

			gsl_blas_dsyr(CblasUpper, gsl_vector_get(model_taumoment, i), &facview1.vector, model_xsigma[i]);
			upperComplete(model_xsigma[i]);
		}
			
		gsl_matrix_memcpy(temp_mat3, model_xsigma[i]);	
		gsl_linalg_LU_decomp(temp_mat3, temp_perm3, &interm1);
		gsl_linalg_LU_invert(temp_mat3, temp_perm3, model_xsigma[i]);

		
		// MEAN UPDATE
		for(j=0; j<n; j++) {
			indview1 = gsl_matrix_column(&dataview.matrix, j);
			gsl_vector_memcpy(temp_vec1, &indview1.vector);
			indview2 = gsl_matrix_row(model_mumean, i);
			gsl_blas_daxpy(-1.0, &indview2.vector, temp_vec1);
			gsl_vector_set_zero(temp_vec2);
			
			facview1 = gsl_matrix_row(model_xmean[i], j);
			gsl_blas_dgemv(CblasTrans, 1.0, model_wmean[i], temp_vec1, 0.0, temp_vec2);

			gsl_blas_dgemv(CblasNoTrans, gsl_vector_get(model_taumoment, i), model_xsigma[i], temp_vec2, 0.0, &facview1.vector);
		}
		
				
	}	
}



void updateZ() {
	// UPDATE Z and associated statistics
	// upperComplete coute cher...
	// + pas de copies de matrices intempestives
	
	for(i=0; i<K; i++) {
		// elements calc for all n
		gsl_matrix_memcpy(temp_mat3, model_wsigma[i]);
		gsl_matrix_scale(temp_mat3, (double)d);
		gsl_blas_dsyrk(CblasUpper, CblasTrans, 1.0, model_wmean[i], 1.0, temp_mat3);
		upperComplete(temp_mat3);
		
		val1 = gsl_sf_psi(gsl_vector_get(model_alpha, i));
		val2 = GSLvectorSum(model_alpha);
		val2 = gsl_sf_psi(val2);
		val1 -= val2;
		

		// ADDITIONAL TERM FOR X
		gsl_matrix_memcpy(temp_mat4, model_xsigma[i]);
		facview1 = gsl_matrix_diagonal(model_xsigma[i]);
		val1 -= 0.5 * GSLvectorSum(&facview1.vector);
		
		val2 = GSLvectorSum(model_musigma[i]);
		indview1 = gsl_matrix_row(model_mumean, i);
		gsl_blas_ddot(&indview1.vector, &indview1.vector, &val3);
		val2 += val3;
		

		val1 = val1 - 0.5 * gsl_vector_get(model_taumoment, i) * val2;


		val4 = gsl_vector_get(model_taumoment, i);
		
		for(j=0; j<n; j++) {

			facview1 = gsl_matrix_row(model_xmean[i], j);
			
			// ADDITONAL TERM FOR X
			gsl_blas_ddot(&facview1.vector, &facview1.vector, &val5);
			val5 *= 0.5;			

			// trace term
			gsl_blas_dsyr(CblasUpper, 1.0, &facview1.vector, temp_mat4);

			gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, temp_mat4, temp_mat3, 0.0, temp_mat5);

			facview2 = gsl_matrix_diagonal(temp_mat5);
			val2 = GSLvectorSum(&facview2.vector);
			// restore temp_mat4
			gsl_blas_dsyr(CblasUpper, -1.0, &facview1.vector, temp_mat4);
			
			
			// cross term
			facview1 = gsl_matrix_row(model_xmean[i], j);
			gsl_blas_dgemv(CblasNoTrans, 1.0, model_wmean[i], &facview1.vector, 0.0, temp_vec1);
			indview1 = gsl_matrix_column(&dataview.matrix, j);
			// decompose for no gsl_vector_memcpy
			gsl_blas_ddot(temp_vec1, &indview1.vector, &val3);
			val2 -= 2.0 * val3;
			
			indview1 = gsl_matrix_row(model_mumean, i);

			gsl_blas_ddot(temp_vec1, &indview1.vector, &val3);
			val2 += 2.0 * val3;
			
			indview1 = gsl_matrix_column(&dataview.matrix, j);
			indview2 = gsl_matrix_row(model_mumean, i);
			gsl_blas_ddot(&indview1.vector, &indview2.vector, &val3);
			val2 = val2 - 2.0 * val3;
			
			val2 += gsl_vector_get(stat_norms, j);
			
			gsl_matrix_set(stat_z, j, i, val1 - val5 - 0.5 * val4 * val2);
		}
	}

	
	// when finished, normalize values for Z
	for(i=0; i<n; i++) {
		grpview1 = gsl_matrix_row(stat_z, i);
		val1 = gsl_vector_max(&grpview1.vector);
		gsl_vector_add_constant(&grpview1.vector, -val1);
		for(j=0; j<K; j++) {
			val1 = gsl_matrix_get(stat_z, i, j);
			val2 = (val1<-700.0) ? 0.0 : gsl_sf_exp(val1);
			gsl_matrix_set(stat_z, i, j, val2);
		}
		
		val1 = GSLvectorSum(&grpview1.vector);
		gsl_vector_scale(&grpview1.vector, 1.0/val1);
	}

}



void updateAlpha() {

	// UPDATE ALPHA
	for(i=0; i<K; i++) {
		gsl_vector_set(model_alpha, i, param_alpha + gsl_vector_get(stat_nk, i));
	}	
}





void updateLB() {
	// terme I
	val1 = 0.0;
	val3 = GSLvectorSum(model_alpha);
	val2 = 0.0;
	for(i=0; i<K; i++) {
		val4 = gsl_vector_get(model_alpha, i);
		val2 += gsl_sf_lngamma(val4);
		val2 += (param_alpha - val4) * (gsl_sf_psi(val4) - gsl_sf_psi(val3));
	}
	val1 += val2 - gsl_sf_lngamma(val3);
	gsl_vector_set(stat_lb, 0, val1);
	
	// terme II
	val1 = 0.0;
	for(i=0; i<K; i++) {
		val2 = gsl_vector_get(model_nua, i);
		for(j=0; j<q; j++) {
			val3 = gsl_matrix_get(model_nub, i, j);
			val1 += (param_nua - val2) * gsl_sf_psi(val2);
			val1 -= param_nua * gsl_sf_log(val3);
			val1 += (1.0 - (param_nub / val3)) * val2;
			val1 += gsl_sf_lngamma(val2);
		}
	}
	gsl_vector_set(stat_lb, 1, val1);
	
	// terme III
	val1 = 0.0;
	for(i=0; i<K; i++) {
		val2 = gsl_vector_get(model_nua, i);
		for(j=0; j<q; j++) {
			val3 = gsl_matrix_get(model_nub, i, j);
			val1 += ((double)d / 2.0) * (gsl_sf_psi(val2) - gsl_sf_log(val3));
		}
		
		gsl_matrix_memcpy(temp_mat3, model_wsigma[i]);
		gsl_matrix_scale(temp_mat3, (double)d);
		gsl_blas_dsyrk(CblasUpper, CblasTrans, 1.0, model_wmean[i], 1.0, temp_mat3);
		upperComplete(temp_mat3);
		gsl_matrix_set_zero(temp_mat4);
		facview1 = gsl_matrix_diagonal(temp_mat4);
		facview2 = gsl_matrix_row(model_numoment, i);
		gsl_vector_memcpy(&facview1.vector, &facview2.vector);
		gsl_blas_dsymm(CblasLeft, CblasUpper, -0.5, temp_mat3, temp_mat4, 0.0, temp_mat5);
		upperComplete(temp_mat5);
		facview1 = gsl_matrix_diagonal(temp_mat5);
		val1 += GSLvectorSum(&facview1.vector);
		
		gsl_matrix_memcpy(temp_mat3, model_wsigma[i]);
		gsl_linalg_LU_decomp(temp_mat3, temp_perm3, &interm1);
		val1 += ((double)d / 2.0) * gsl_linalg_LU_lndet(temp_mat3);
	}
	gsl_vector_set(stat_lb, 2, val1);
	
	// terme IV
	val1 = 0.0;
	for(i=0; i<K; i++) {
		val2 = GSLvectorSum(model_musigma[i]);
		indview1 = gsl_matrix_row(model_mumean, i);
		gsl_blas_ddot(&indview1.vector, &indview1.vector, &val3);
		val2 += val3;
		val2 *= -param_nustar / 2.0;
		val1 += val2;
		
		indview2 = gsl_matrix_row(param_mustar, i);
		gsl_blas_ddot(&indview1.vector, &indview2.vector, &val2);
		val1 += param_nustar * val2;
		
		for(j=0; j<d; j++) {

			val1 += 0.5 * gsl_sf_log(gsl_vector_get(model_musigma[i], j));
		}
	}
	gsl_vector_set(stat_lb, 3, val1);
	
	// terme V
	val1 = 0.0;
	for(i=0; i<K; i++) {
		for(j=0; j<n; j++) {
			val2 = gsl_matrix_get(stat_z, j, i);
			if(val2 > 0.0) {
				val1 -= val2 * gsl_sf_log(val2);
			}
		}
		
		val2 = gsl_vector_get(model_alpha, i);
		val3 = GSLvectorSum(model_alpha);
		val1 += gsl_vector_get(stat_nk, i) * (gsl_sf_psi(val2) - gsl_sf_psi(val3));
	}
	gsl_vector_set(stat_lb, 4, val1);
	
	// terme VI
	val1 = 0.0;
	for(i=0; i<K; i++) {
		facview1 = gsl_matrix_diagonal(stat_Sk[i]);
		val1 -= 0.5 * GSLvectorSum(&facview1.vector);
		
		gsl_matrix_memcpy(temp_mat3, model_xsigma[i]);
		gsl_linalg_LU_decomp(temp_mat3, temp_perm3, &interm1);
		val1 += 0.5 * gsl_vector_get(stat_nk, i) * gsl_linalg_LU_lndet(temp_mat3);
	}
	gsl_vector_set(stat_lb, 5, val1);
	
	// terme VII
	val1 = 0.0;
	for(i=0; i<K; i++) {
		val2 = gsl_vector_get(stat_yk2, i);
		
		gsl_matrix_memcpy(temp_mat3, model_wsigma[i]);
		gsl_matrix_scale(temp_mat3, (double)d);
		gsl_blas_dsyrk(CblasUpper, CblasTrans, 1.0, model_wmean[i], 1.0, temp_mat3);
		upperComplete(temp_mat3);
		gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, temp_mat3, stat_Sk[i], 0.0, temp_mat5);
		upperComplete(temp_mat5);
		facview1 = gsl_matrix_diagonal(temp_mat5);
		val2 += GSLvectorSum(&facview1.vector);
		
		val3 = GSLvectorSum(model_musigma[i]);
		indview1 = gsl_matrix_row(model_mumean, i);
		gsl_blas_ddot(&indview1.vector, &indview1.vector, &val4);
		val3 += val4;
		val2 += gsl_vector_get(stat_nk, i) * val3;
		
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, stat_syk[i], model_wmean[i], 0.0, temp_mat3);
		facview1 = gsl_matrix_diagonal(temp_mat3);
		val2 -= 2.0 * GSLvectorSum(&facview1.vector);
		
		facview1 = gsl_matrix_row(stat_sk, i);
		gsl_blas_dgemv(CblasNoTrans, 1.0, model_wmean[i], &facview1.vector, 0.0, temp_vec1);
		gsl_blas_ddot(&indview1.vector, temp_vec1, &val3);
		val2 += 2.0 * val3;
		
		indview2 = gsl_matrix_row(stat_yk, i);
		gsl_blas_ddot(&indview2.vector, &indview1.vector, &val3);
		val2 -= 2.0 * val3;


		val1 -= 0.5 * gsl_vector_get(model_taumoment, i) * val2;
	}
	gsl_vector_set(stat_lb, 6, val1);


	
	
	// terme VIII
	val1 = 0.0;
	/*
	for(i=0; i<K; i++) {
		
		val2 = gsl_vector_get(model_taua, i);
		val3 = gsl_vector_get(model_taub, i);
		val4 = gsl_vector_get(model_taumoment, i);
		
		val1 += gsl_sf_lngamma(val2);
		val1 += (param_taua - val2) * gsl_sf_psi(val2);
		val1 -= param_taua * gsl_sf_log(val3);
		val1 += val2 - param_taub * val4;
	}
	*/
	gsl_vector_set(stat_lb, 7, val1);
}




void updateMu() {
	for(i=0; i<K; i++) {
		// COVARIANCE UPDATE
		gsl_vector_set_all(model_musigma[i], param_nustar + gsl_vector_get(model_taumoment, i) * gsl_vector_get(stat_nk, i));
		
		for(j=0; j<d; j++) {
			val1 = gsl_vector_get(model_musigma[i], j);
			val1 = 1.0/val1;
			gsl_vector_set(model_musigma[i], j, val1);
		}
		
		// MEAN UPDATE
		gsl_vector_set_zero(temp_vec1);

		indview1 = gsl_matrix_row(stat_yk, i);
		gsl_blas_daxpy(1.0, &indview1.vector, temp_vec1);
		facview1 = gsl_matrix_row(stat_sk, i);
		gsl_blas_dgemv(CblasNoTrans, -1.0, model_wmean[i], &facview1.vector, 1.0, temp_vec1);
		
		//
		gsl_vector_scale(temp_vec1, gsl_vector_get(model_taumoment, i));
		
		
		indview1 = gsl_matrix_row(param_mustar, i);
		gsl_blas_daxpy(param_nustar, &indview1.vector, temp_vec1);
		indview1 = gsl_matrix_row(model_mumean, i);
		gsl_vector_mul(temp_vec1, model_musigma[i]);
		gsl_vector_memcpy(&indview1.vector, temp_vec1);
	}
}



void updateW() {
	for(i=0; i<K; i++) {
		// COVARIANCE UPDATE
		gsl_matrix_set_zero(temp_mat4);
		facview1 = gsl_matrix_row(model_numoment, i);
		facview2 = gsl_matrix_diagonal(temp_mat4);
		gsl_vector_memcpy(&facview2.vector, &facview1.vector);
		
		gsl_matrix_memcpy(temp_mat5, stat_Sk[i]);
		gsl_matrix_scale(temp_mat5, gsl_vector_get(model_taumoment, i));
		gsl_matrix_add(temp_mat4, temp_mat5);
		
		gsl_linalg_LU_decomp(temp_mat4, temp_perm3, &interm1);
		gsl_linalg_LU_invert(temp_mat4, temp_perm3, model_wsigma[i]);
		
		// MEAN UPDATE
		indview1 = gsl_matrix_row(model_mumean, i);
		for(j=0; j<d; j++) {
			facview1 = gsl_matrix_row(stat_syk[i], j);
			gsl_vector_memcpy(temp_vec2, &facview1.vector);
			facview1 = gsl_matrix_row(stat_sk, i);
			gsl_blas_daxpy(-gsl_vector_get(&indview1.vector, j), &facview1.vector, temp_vec2);

			facview1 = gsl_matrix_row(model_wmean[i], j);
			gsl_blas_dgemv(CblasNoTrans, gsl_vector_get(model_taumoment, i), model_wsigma[i], temp_vec2, 0.0, &facview1.vector);
		}
	}
}



void updateNu() {
	// UPDATE NU
	for(i=0; i<K; i++) {
		gsl_vector_set(model_nua, i, param_nua + ((double)d/2.0));
		

		for(j=0; j<d; j++) {
			gsl_matrix_memcpy(temp_mat3, model_wsigma[i]);
			facview1 = gsl_matrix_row(model_wmean[i], j);
			gsl_blas_dsyr(CblasUpper, 1.0, &facview1.vector, temp_mat3);
			upperComplete(temp_mat3);
			facview1 = gsl_matrix_diagonal(temp_mat3);
			facview2 = gsl_matrix_row(stat_trace, j);
			gsl_vector_memcpy(&facview2.vector, &facview1.vector);
		}

		facview1 = gsl_matrix_row(model_nub, i);
		facview2 = gsl_matrix_row(model_numoment, i);
		val2 = gsl_vector_get(model_nua, i);
		
		for(j=0; j<q; j++) {
			indview1 = gsl_matrix_column(stat_trace, j);
			val1 = GSLvectorSum(&indview1.vector);
			val1 = param_nub + 0.5 * val1;
			gsl_vector_set(&facview1.vector, j, val1);
			gsl_vector_set(&facview2.vector, j, val2 / val1);
		}
	}	
}




void updateTau() {
	for(i=0; i<K; i++) {
		val2 = gsl_vector_get(stat_yk2, i);
		
		gsl_matrix_memcpy(temp_mat3, model_wsigma[i]);
		gsl_matrix_scale(temp_mat3, (double)d);
		gsl_blas_dsyrk(CblasUpper, CblasTrans, 1.0, model_wmean[i], 1.0, temp_mat3);
		upperComplete(temp_mat3);
		gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, temp_mat3, stat_Sk[i], 0.0, temp_mat5);
		upperComplete(temp_mat5);
		facview1 = gsl_matrix_diagonal(temp_mat5);
		val2 += GSLvectorSum(&facview1.vector);
		
		val3 = GSLvectorSum(model_musigma[i]);
		indview1 = gsl_matrix_row(model_mumean, i);
		gsl_blas_ddot(&indview1.vector, &indview1.vector, &val4);
		val3 += val4;
		val2 += gsl_vector_get(stat_nk, i) * val3;
		
		
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, stat_syk[i], model_wmean[i], 0.0, temp_mat3);
		facview1 = gsl_matrix_diagonal(temp_mat3);
		val2 -= 2.0 * GSLvectorSum(&facview1.vector);
		
		facview1 = gsl_matrix_row(stat_sk, i);
		gsl_blas_dgemv(CblasNoTrans, 1.0, model_wmean[i], &facview1.vector, 0.0, temp_vec1);
		gsl_blas_ddot(&indview1.vector, temp_vec1, &val3);
		val2 += 2.0 * val3;
		
		indview2 = gsl_matrix_row(stat_yk, i);
		gsl_blas_ddot(&indview2.vector, &indview1.vector, &val3);
		val2 -= 2.0 * val3;

		gsl_vector_set(model_taua, i, param_taua + gsl_vector_get(stat_nk, i) * (double)d * 0.5);
		gsl_vector_set(model_taub, i, param_taub + val2 * 0.5);
		gsl_vector_set(model_taumoment, i, gsl_vector_get(model_taua, i) / gsl_vector_get(model_taub, i));
	}		
}



void displayCurModel() {
	interm1 = 0;
	for(i=0; i<K; i++) {
		val2 = gsl_vector_get(model_alpha, i);
		if(val2 > 0.002) {
			interm1++;
		}
	}
	Rprintf("%d active comp.\n", interm1);
}



void displayBound(int ind) {
	val1 = GSLvectorSum(stat_lb);
	if(!ind) {
		Rprintf("bound is %f ", val1);
	} else {
		Rprintf("it %d, bound is %f ", ind, val1);
	}
	
	/*
	for(i=0; i<stat_lb->size; i++) {
		Rprintf("%f ", gsl_vector_get(stat_lb, i));
	}
	Rprintf("\n");
	*/
	
}


void updateSuffStats() {
	
	for(i=0; i<K; i++) {
		// Nk & yk2
		varview1 = gsl_matrix_column(stat_z, i);
		val1 = GSLvectorSum(&varview1.vector);
		gsl_vector_set(stat_nk, i, val1);
		gsl_blas_ddot(&varview1.vector, stat_norms, &val2);
		gsl_vector_set(stat_yk2, i, val2);
		
		// yk & sk & syk & Sk
		indview1 = gsl_matrix_row(stat_yk, i);
		facview1 = gsl_matrix_row(stat_sk, i);
		gsl_vector_set_zero(&indview1.vector);
		gsl_vector_set_zero(&facview1.vector);
		gsl_matrix_set_zero(stat_syk[i]);
		gsl_matrix_memcpy(stat_Sk[i], model_xsigma[i]);
		gsl_matrix_scale(stat_Sk[i], val1);
		
		gsl_vector *suffStatsTemp;
		suffStatsTemp = &varview1.vector;
		
		for(j=0; j<n; j++) {

			val2 = suffStatsTemp->data[j * suffStatsTemp->stride];
			indview2 = gsl_matrix_column(&dataview.matrix, j);
			gsl_blas_daxpy(val2, &indview2.vector, &indview1.vector);
			facview2 = gsl_matrix_row(model_xmean[i], j);
			gsl_blas_daxpy(val2, &facview2.vector, &facview1.vector);
			gsl_blas_dger(val2, &indview2.vector, &facview2.vector, stat_syk[i]);
			gsl_blas_dsyr(CblasUpper, val2, &facview2.vector, stat_Sk[i]);

		}
		upperComplete(stat_Sk[i]);
	}
}


// method = "sample" or "choice"
// - "sample" : classic init
// - "choice" : init means with randomly chosen elements (n must be > to ncomp)
// - "partial" : output model for further input. (default with init method "sample" ; should be allowed to be given a numeric vector in future)
//SEXP mppca(SEXP data, SEXP ncomp, SEXP maxit, SEXP qmax, SEXP method, SEXP inmod, SEXP inw) {
SEXP mppca(SEXP data, SEXP ncomp, SEXP thres, SEXP maxit, SEXP qmax) {
	
	Rprintf("allocating...\n");
	alloc(data, ncomp, qmax);
	
	Rprintf("initializing...\n");
	init();
	
	/*
	if(!strcmp(CHAR(STRING_ELT(method, 0)), "choice")) {
		initWithElements();
	}
	*/
	
	/*
	if(inmod != R_NilValue) {
		initWithReadyModel(inmod);
	}
	*/
	
	/*
	if(inw != R_NilValue) {
		initWithReadyW(inw);
	}
	*/

	Rprintf("starting...\n");

	updateZ();
	updateSuffStats();
	
	updateMu();
	updateX();
	updateSuffStats();


	updateLB();
	// compute LB
	lb = GSLvectorSum(stat_lb);
	lbvar = GSL_POSINF;
	nit = 1;
	
	
	SEXP stamp, newstamp;
	
	double *local_bounds = 0;
	unsigned converge = 0;
	
	while(!converge) {
		PROTECT(stamp=getTimestamp());
		displayBound(nit);
		
		if(local_bounds==0) {
			local_bounds = Calloc(1, double);
			local_bounds[0] = lb;
		} else {
			local_bounds = Realloc(local_bounds, nit, double);
			local_bounds[nit-1] = lb;
		}
		
		displayCurModel();

		updateAlpha();
		
		for(dummy=0; dummy<1; dummy++) {
			updateX();
			updateZ();
		}
		updateSuffStats();
		
				


		//Rprintf("pass");
		for(dummy=0; dummy<1; dummy++) {
		//	Rprintf(" %d", dummy+1);
			updateMu();
			updateW();
			updateNu();
			
		}
		//Rprintf("\n");
		
		

		updateLB();		
		// COMPUTE LB
		newlb = GSLvectorSum(stat_lb);
		lbvar = newlb - lb;
		lb = newlb;
		nit++;

		if(maxit==R_NilValue) {
			if(lbvar < REAL(coerceVector(thres, REALSXP))[0]) converge = 1;
		} else {
			if(nit > INTEGER(coerceVector(maxit, INTSXP))[0]) converge = 1;
		}

		PROTECT(newstamp=getElapsed(stamp));
		//Rprintf("step lasted %fs\n", REAL(newstamp)[0]);
		UNPROTECT(2);

		
	}

	Rprintf("final bound was %f\n", lb);
	/*
	for(i=0; i<stat_lb->size; i++) {
		Rprintf("%f ", gsl_vector_get(stat_lb, i));
	}
	Rprintf("\n");
	*/

	if(local_bounds==0) {
		local_bounds = Calloc(1, double);
		local_bounds[0] = lb;
	} else {
		local_bounds = Realloc(local_bounds, nit, double);
		local_bounds[nit-1] = lb;
	}

//	if(strcmp(CHAR(STRING_ELT(method, 0)), "partial")) {
	
	
		Rprintf("eigen mppca...\n");
	
		// After Main loop, perform eigen MPPCA
		for(i=0; i<K; i++) {
			val1 = gsl_vector_get(model_alpha, i);
			if(val1 > 2.0) {
				// try to 0 all values below some thres (say 10^(-30))
				
				gsl_blas_dsyrk(CblasUpper, CblasTrans, 1.0, model_wmean[i], 0.0, temp_mat3);
				upperComplete(temp_mat3);
				
				for(j=0; j<q; j++) {
					for(k=0; k<q; k++) {
						val1 = gsl_matrix_get(temp_mat3, j, k);
						if(fabs(val1) < pow(10.0, -30.0)) gsl_matrix_set(temp_mat3, j, k, 0.0);
					}
				}
				
				gsl_eigen_symmv(temp_mat3, temp_vec2, temp_mat4, temp_symmv1);

				// eigenvalues and vectors are unordered : must be reordered.
		
				gsl_sort_vector_index(temp_perm3, temp_vec2);
				
				// permutation : index of vector elt which would have been stored here, in asc order.
				for(j=0; j<q; j++) {
					facview1 = gsl_matrix_column(temp_mat3, j);
					facview2 = gsl_matrix_column(temp_mat4, temp_perm3->data[q-1-j]);
					gsl_vector_memcpy(&facview1.vector, &facview2.vector);
				}
				// when done, post multiply
				gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, model_wmean[i], temp_mat3, 0.0, temp_mat6);
				gsl_matrix_memcpy(model_wmean[i], temp_mat6);
			}
		}
	
		Rprintf("post-processing...\n");
		
		updateNu();
	
		for(dummy=0; dummy<1; dummy++) {
			//Rprintf("pass %d...\n", dummy+1);
		
			updateX();
		}
	
		
		// before building R structures, make all factors the same sign
		for(i=0; i<K; i++) {
			for(j=0; j<q; j++) {
				indview1 = gsl_matrix_column(model_wmean[i], j);
				if(gsl_vector_get(&indview1.vector, 0) < 0.0) {
					gsl_blas_dscal(-1.0, &indview1.vector);
				}
			}
		}
		

//	}



	// allocate and build returned data
	SEXP model_out, alpha_out, numoment_out, nua_out, nub_out, taumoment_out, taua_out, taub_out, wmean_out, wsigma_out, xmean_out, xsigma_out, musigma_out, mumean_out;
	SEXP mustar_out;
	SEXP curnumoment, curnub, curwmean, curwsigma, curxmean, curxsigma, curmumean, curmusigma, curmustar;
	SEXP dims1;

	SEXP Rbound;
	PROTECT(alpha_out = allocVector(REALSXP, K));
	PROTECT(numoment_out = allocVector(VECSXP, K));
	PROTECT(nua_out = allocVector(REALSXP, K));
	PROTECT(nub_out = allocVector(VECSXP, K));
	PROTECT(taumoment_out = allocVector(REALSXP, K));
	PROTECT(taua_out = allocVector(REALSXP, K));	
	PROTECT(taub_out = allocVector(REALSXP, K));
	PROTECT(wmean_out = allocVector(VECSXP, K));
	PROTECT(wsigma_out = allocVector(VECSXP, K));
	PROTECT(xmean_out = allocVector(VECSXP, K));
	PROTECT(xsigma_out = allocVector(VECSXP, K));
	PROTECT(mumean_out = allocVector(VECSXP, K));
	PROTECT(mustar_out = allocVector(VECSXP, K));
	PROTECT(musigma_out = allocVector(VECSXP, K));
	//
	PROTECT(Rbound = allocVector(REALSXP, nit));

	
	for(i=0; i<K; i++) {
		SET_VECTOR_ELT(numoment_out, i, curnumoment = allocVector(REALSXP, q));
		SET_VECTOR_ELT(nub_out, i, curnub = allocVector(REALSXP, q));
		SET_VECTOR_ELT(wmean_out, i, curwmean = allocMatrix(REALSXP, d, q));
		SET_VECTOR_ELT(wsigma_out, i, curwsigma = allocMatrix(REALSXP, q,q));
		SET_VECTOR_ELT(xmean_out, i, curxmean = allocMatrix(REALSXP, n, q));
		SET_VECTOR_ELT(xsigma_out, i, curxsigma = allocMatrix(REALSXP, q,q));
		SET_VECTOR_ELT(mumean_out, i, curmumean = allocVector(REALSXP, d));
		//SET_VECTOR_ELT(musigma_out, i, curmusigma = allocMatrix(REALSXP, d,d));
		SET_VECTOR_ELT(musigma_out, i, curmusigma = allocVector(REALSXP,d));
		SET_VECTOR_ELT(mustar_out, i, curmustar = allocVector(REALSXP, d));

		
		REAL(alpha_out)[i] = gsl_vector_get(model_alpha, i);
		REAL(nua_out)[i] = gsl_vector_get(model_nua, i);
		REAL(taumoment_out)[i] = gsl_vector_get(model_taumoment, i);
		REAL(taua_out)[i] = gsl_vector_get(model_taua, i);
		REAL(taub_out)[i] = gsl_vector_get(model_taub, i);


		for(j=0; j<d; j++) {
			REAL(curmumean)[j] = gsl_matrix_get(model_mumean, i, j);
			REAL(curmustar)[j] = gsl_matrix_get(param_mustar, i, j);
			/*
			for(k=0; k<d; k++) {
				REAL(curmusigma)[k*d + j] = gsl_matrix_get(model_musigma[i], j, k);
			}*/
			REAL(curmusigma)[j] = gsl_vector_get(model_musigma[i], j);
			for(k=0; k<q; k++) {
				REAL(curwmean)[k*d + j] = gsl_matrix_get(model_wmean[i], j, k);
			}
		}
		

		for(j=0; j<q; j++) {
			REAL(curnub)[j] = gsl_matrix_get(model_nub, i, j);
			REAL(curnumoment)[j] = gsl_matrix_get(model_numoment, i, j);
			for(k=0; k<q; k++) {
				REAL(curwsigma)[k*q + j] = gsl_matrix_get(model_wsigma[i], j, k);
				REAL(curxsigma)[k*q + j] = gsl_matrix_get(model_xsigma[i], j, k);
			}
			
			for(k=0; k<n; k++) {
				REAL(curxmean)[j*n + k] = gsl_matrix_get(model_xmean[i], k, j);
			}
			
		}
		
	}
	
	for(i=0; i<nit; i++) {
		REAL(Rbound)[i] = local_bounds[i];
	}

	Free(local_bounds);
	
	PROTECT(model_out = allocVector(VECSXP, 14));
	SET_VECTOR_ELT(model_out, 0, alpha_out);
	SET_VECTOR_ELT(model_out, 1, numoment_out);
	SET_VECTOR_ELT(model_out, 2, nua_out);
	SET_VECTOR_ELT(model_out, 3, nub_out);	
	SET_VECTOR_ELT(model_out, 4, taumoment_out);
	SET_VECTOR_ELT(model_out, 5, taua_out);
	SET_VECTOR_ELT(model_out, 6, taub_out);
	SET_VECTOR_ELT(model_out, 7, wmean_out);
	SET_VECTOR_ELT(model_out, 8, wsigma_out);
	SET_VECTOR_ELT(model_out, 9, xmean_out);
	SET_VECTOR_ELT(model_out, 10, xsigma_out);
	SET_VECTOR_ELT(model_out, 11, mumean_out);
	SET_VECTOR_ELT(model_out, 12, musigma_out);
	SET_VECTOR_ELT(model_out, 13, mustar_out);
	SET_VECTOR_ELT(model_out, 14, Rbound);

	
	
	PROTECT(dims1 = allocVector(VECSXP, 14));
	SET_VECTOR_ELT(dims1, 0, mkChar("alpha"));
	SET_VECTOR_ELT(dims1, 1, mkChar("numoment"));
	SET_VECTOR_ELT(dims1, 2, mkChar("nua"));
	SET_VECTOR_ELT(dims1, 3, mkChar("nub"));
	SET_VECTOR_ELT(dims1, 4, mkChar("taumoment"));
	SET_VECTOR_ELT(dims1, 5, mkChar("taua"));
	SET_VECTOR_ELT(dims1, 6, mkChar("taub"));
	SET_VECTOR_ELT(dims1, 7, mkChar("wmean"));
	SET_VECTOR_ELT(dims1, 8, mkChar("wsigma"));
	SET_VECTOR_ELT(dims1, 9, mkChar("xmean"));	
	SET_VECTOR_ELT(dims1, 10, mkChar("xsigma"));
	SET_VECTOR_ELT(dims1, 11, mkChar("mumean"));
	SET_VECTOR_ELT(dims1, 12, mkChar("musigma"));
	SET_VECTOR_ELT(dims1, 13, mkChar("mustar"));
	SET_VECTOR_ELT(dims1, 14, mkChar("bound"));
	setAttrib(model_out, R_NamesSymbol, dims1);


	end();
	UNPROTECT(17);	

	return(model_out);
		
}


SEXP getResp(SEXP data, SEXP model) {
	K = length(getListElement(model, "alpha"));
	q = length(VECTOR_ELT(getListElement(model, "numoment"), 0));
	SEXP R_K, R_q;
	PROTECT(R_K=allocVector(INTSXP, 1));
	PROTECT(R_q=allocVector(INTSXP, 1));
	INTEGER(R_K)[0] = K;
	INTEGER(R_q)[0] = q;
	
	alloc(data, R_K, R_q);
	init();
	initWithReadyModel(model);
	
	updateZ();
	SEXP z_out;
	PROTECT(z_out = allocMatrix(REALSXP, n, K));
	for(i=0; i<n; i++) {
		for(j=0; j<K; j++) {
			REAL(z_out)[j*n + i] = gsl_matrix_get(stat_z, i, j);
		}
	}
	
	end();
	UNPROTECT(3);
	return(z_out);
}


