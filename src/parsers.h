// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#ifndef PARSERS_H
#define PARSERS_H

#include <R.h>
#include <Rinternals.h>

typedef struct {
	int nelts;
	int *set;
	double *setdists;
} summary;

#ifdef __cplusplus
extern "C" {
#endif

void transformPath(char *);
char *writeUrls(summary, SEXP, char *);

#ifdef __cplusplus
}
#endif

#endif
