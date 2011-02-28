// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#ifndef MMPPCA_H
#define MMPPCA_H

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

// declare all variables
int K,L,d,q,i,j,k, a, b, interm1,nit,dummy, c_maxit;
double val1,val2,val3, val4, val5, val6, val7, val8, val9, val10, lb, newlb, lbvar;
double param_alpha, param_nua, param_nub, param_numoment, param_taua, param_taub, param_taumoment, param_nustar;

SEXP in_alpha, in_wmean, in_mumean;

gsl_vector_view varview1;
gsl_vector_view varview2;
gsl_vector_view facview1;
gsl_vector_view facview2;
gsl_vector_view indview1;
gsl_vector_view indview2;
gsl_vector_view grpview1;
gsl_vector_view grpview2;
gsl_matrix_view tempview1;


gsl_matrix *in_smp;
gsl_matrix **in_smpfactors;
gsl_vector *in_nonvoidfactors;
gsl_vector *in_scales1;
gsl_vector *in_scales2;

gsl_matrix *stat_z;
gsl_vector *stat_norms1;
gsl_vector *stat_norms2;
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
gsl_matrix **model_x1mean;
gsl_matrix ***model_x2mean;
gsl_matrix **model_wsigma;

gsl_matrix *temp_mat1;
gsl_matrix *temp_mat2;
gsl_matrix *temp_mat3;
gsl_matrix *temp_mat4;
gsl_matrix *temp_mat5;
gsl_matrix *temp_mat6;
gsl_permutation *temp_perm1;
gsl_permutation *temp_perm3;
gsl_permutation *temp_perm4;
gsl_vector *temp_vec1;
gsl_vector *temp_vec2;
gsl_vector *temp_vec3;
gsl_vector *temp_vec4;
gsl_vector *temp_vec5;
gsl_vector *temp_vec6;
gsl_eigen_symmv_workspace *temp_symmv1;

double *mins;
double *maxs;


// functions

int isNonVoid(gsl_vector *);

void allocM(SEXP, SEXP);
void initM();
void initMWithPrior(SEXP, SEXP);
void initMWithReadyModel(SEXP);
void endM();

SEXP mmppca(SEXP, SEXP, SEXP, SEXP);

void updateMX();
void updateMZ();
void updateMLB();
void updateMMu();
void updateMW();
void updateMAlpha();
void updateMNu();
void updateMSuffStats();
void displayMCurModel();
void printDiagnostic();
void setXsampleCov();
void updateMTau();




#endif
