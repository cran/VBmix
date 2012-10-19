// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#ifndef UTILS_H
#define UTILS_H




#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <sys/time.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

void printSxpIntVector(SEXP);
void printSxpRealVector(SEXP);
void printSxpMatrix(SEXP);
void printGSLvector(gsl_vector *);
void printGSLmatrix(gsl_matrix *);
void printDoubleVector(double *, int);

int whichmax(gsl_vector *);
SEXP getListElement(SEXP, const char *);
void vectorToSXP(SEXP *, gsl_vector *);
void intVectorToSXP(SEXP *, gsl_vector *);
void matrixToSXP(SEXP *, gsl_matrix *);
void SXPtoVector(gsl_vector *, SEXP);
void SXPtoMatrix(gsl_matrix *, SEXP);
double SXPvectorSum(SEXP);
double GSLvectorSum(gsl_vector *);
SEXP multinomial(SEXP, SEXP);
SEXP mvngen(SEXP, SEXP, SEXP);
SEXP gmmgen(SEXP, SEXP);
SEXP mvndensity(SEXP, SEXP, SEXP);
SEXP gmmdensity(SEXP, SEXP);
SEXP klmc(SEXP, SEXP, SEXP);
SEXP jsmc(SEXP, SEXP, SEXP);
SEXP klut(SEXP, SEXP);
SEXP jsut(SEXP, SEXP);
SEXP extractSimpleModel(SEXP, SEXP);
SEXP buildPlainMod(SEXP, SEXP, SEXP);
SEXP appendRealVector(SEXP, SEXP);
SEXP appendList(SEXP, SEXP);

SEXP GSLklut(SEXP, SEXP);
SEXP GSLklutSet(SEXP, SEXP);


SEXP getTimestamp();
SEXP getElapsed(SEXP);


// double usage lists (wasteful but convenient) => stores int and double at the same time
typedef struct dummy {
	int intVal;
	double doubleVal;
	struct dummy *next;
} item;

typedef struct {
	item *first;
	item *last;
} list;

list *createList();
void addItem(list *, int, double);
void insertItem(list *, int, double, int);
int size(list *);
int getIntItem(list *, int);
void setDoubleItem(list *, int, double);
double getDoubleItem(list *, int);
void dropItem(list *, int);
void dropList(list *);
void printList(list *);

typedef struct metadummy {
	int key;
	list *lst;
	struct metadummy *next;
} metaitem;

typedef struct {
	metaitem *first;
	metaitem *last;
} metalist;

void addList(metalist *, int);
int metaSize(metalist *);
list *getList(metalist *, int);
void dropMetalist(metalist *);
int getIntMetaItem(metalist *, int);


SEXP listToSXP(list *);


// test presence of a specific int value
// and index of it
int indexOfInt(list *, int);
// searches int items in first elements of lists in the metalist
int indexOfKey(metalist *, int);


void listToIntVector(gsl_vector *, list *);
void listToDoubleVector(gsl_vector *, list *);

double sum(double *, int, int);


SEXP pointwise(SEXP, SEXP);
int gramschmidt(gsl_matrix *);
SEXP R_gramschmidt(const SEXP);
void upperComplete(gsl_matrix *);
SEXP sort_index(SEXP, SEXP);

void GSLsample(int, int, int*);
SEXP sample(SEXP, SEXP);
int contains(int, int, int*);

SEXP Rdct(SEXP);
SEXP Rdct2D(SEXP);
SEXP RinvDct2D(SEXP);
SEXP control(SEXP, SEXP);
double getDistance(gsl_vector *, gsl_vector *);
void getColumnNorms(gsl_matrix *, gsl_vector *);
void getCovariance(gsl_matrix *, gsl_matrix *);
void getMean(gsl_matrix *, gsl_vector *);

SEXP rDirichlet(SEXP, SEXP);

#ifdef __cplusplus
}
#endif

#endif





