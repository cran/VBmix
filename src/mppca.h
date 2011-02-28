// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#ifndef MPPCA_H
#define MPPCA_H

#include <QtCore>
#include <iostream>

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
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "utils.h"

// declare all variables
int K,n,d,q,i,j,k,interm1,nit,dummy, lb_store;
double val1,val2,val3, val4, val5, val6, lb, newlb, lbvar;
double param_alpha, param_nua, param_nub, param_numoment, param_taua, param_taub, param_taumoment, param_nustar;
gsl_vector_view varview1;
gsl_vector_view varview2;
gsl_vector_view facview1;
gsl_vector_view facview2;
gsl_vector_view indview1;
gsl_vector_view indview2;
gsl_vector_view grpview1;
gsl_vector_view grpview2;
gsl_matrix_view dataview;
gsl_matrix_view tempview1;
gsl_matrix_view tempview2;

gsl_matrix *stat_z;
gsl_vector *stat_norms;
gsl_vector *stat_lb;
gsl_matrix *stat_trace;
gsl_vector *stat_nk;
gsl_matrix *stat_yk;
gsl_matrix *stat_sk;
gsl_matrix **stat_syk;
gsl_matrix **stat_Sk;
gsl_vector *stat_yk2;

gsl_matrix *param_mustar;
gsl_vector *model_alpha;
gsl_vector *model_nua;
gsl_matrix *model_nub;
gsl_matrix *model_numoment;
gsl_vector *model_taua;
gsl_vector *model_taub;
gsl_vector *model_taumoment;
gsl_matrix *model_mumean;
gsl_vector **model_musigma;
gsl_matrix **model_wmean;
gsl_matrix **model_xsigma;
gsl_matrix **model_xmean;
gsl_matrix **model_wsigma;

gsl_matrix *temp_mat1;
gsl_matrix *temp_mat2;
gsl_matrix *temp_mat3;
gsl_matrix *temp_mat4;
gsl_matrix *temp_mat5;
gsl_matrix *temp_mat6;
gsl_matrix *temp_mat7;
gsl_permutation *temp_perm1;
gsl_permutation *temp_perm3;
gsl_permutation *temp_perm4;
gsl_permutation *temp_perm5;
gsl_vector *temp_vec1;
gsl_vector *temp_vec2;
gsl_vector *temp_vec3;
gsl_vector *temp_vec4;
gsl_vector *temp_vec5;
gsl_eigen_symmv_workspace *temp_symmv1;

double *mins;
double *maxs;

SEXP output;

// functions
extern "C" {
	void alloc(SEXP, SEXP, SEXP);
	void init();
	void initWithPrior(SEXP, SEXP, SEXP);
	//void initWithElements();
	void initWithReadyModel(SEXP);
	void initWithReadyW(SEXP);
	void end();

	SEXP mppca(SEXP, SEXP, SEXP, SEXP, SEXP);

	void updateX();
	void updateZ();
	void updateLB();
	void updateMu();
	void updateW();
	void updateAlpha();
	void updateNu();
	void displayCurModel();
	void displayBound(int);
	void updateSuffStats();

	void writeMatrix(gsl_matrix *);
	
	SEXP getResp(SEXP, SEXP);
}


#endif
