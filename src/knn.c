// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "utils.h"
#include "vbconstr.h"
#include "vbcomp.h"
#include "varbayes.h"
#include <stdlib.h>

typedef struct cell {
	int ind;
	double val;
} cell;

int compare(const void *a, const void *b) {
	if(((cell *)a)->val > ((cell *)b)->val) {
		return 1;
	} else if (((cell *)a)->val < ((cell *)b)->val) {
		return -1;
	} else {
		return 0;
	}
}



SEXP mixKnn(SEXP data, SEXP labels, SEXP n, SEXP KLparam) {
	// runs the knn algorithm on leave-1-out fashion with the specified n
	// use 500 for jsmc
	
	
	int i,j;
	int errcount;
	int len = length(labels);
	int found;
	cell dists[len];
	//list *kcounts = createList();
	list *kcounts = Calloc(1, list);
	item *ptr;
	int cur;
	
	// coerce vectors appropriately
	PROTECT(labels = coerceVector(labels, INTSXP));
	PROTECT(n = coerceVector(n, INTSXP));
	
	for(i=0; i<length(labels); i++) {
		if(indexOfInt(kcounts, INTEGER(labels)[i]) == -1) {
			addItem(kcounts, INTEGER(labels)[i], 0.0);
		}
	}
		
	
	SEXP nits;
	PROTECT(nits = allocVector(INTSXP, 1));
	INTEGER(nits)[0] = 500;
	
	// get levels of the labels vector
	// a priori, the largest value in the vector is the maximum number of levels.
	// check if at least one value of each level is found (max and is.found functions,
	// chained list to store found values)
	
	
	
	for(i=0; i<len; i++) {
		Rprintf("element %d...\n", i);
		// first, compute jsmc wrt to each other element in data
		Rprintf("computing distances...\n");
		for(j=0; j<len; j++) {
			dists[j].ind = j;
			dists[j].val = REAL(jsmc(VECTOR_ELT(data, i), VECTOR_ELT(data, j), KLparam))[0];
		}
		
		// then sort 
		qsort(dists, len, sizeof(cell), compare);
		
		// knn search : search the first label which has n closest elements
		found=0;
		j=1;
		
		Rprintf("finding closest label...\n");
		while(!found) {
			cur = INTEGER(labels)[dists[j].ind];
			
			Rprintf("seek label %d from item %d\n", cur, dists[j].ind);
			
			setDoubleItem(kcounts, indexOfInt(kcounts, cur), getDoubleItem(kcounts, indexOfInt(kcounts, cur)) +1.0);                                
			
			Rprintf("cur: %d, count: %f\n", cur, getDoubleItem(kcounts, indexOfInt(kcounts, cur)));
			
			
			if(getDoubleItem(kcounts, indexOfInt(kcounts, cur)) >= INTEGER(n)[0]) {
				found = 1;
			}
			j++;
			
		}
		
		
		if(cur != INTEGER(labels)[i]) {
			errcount +=1;
			Rprintf("mismatch...\n");
		} else {
			Rprintf("success.\n");
		}
		
		// reset kcount structure to zeroes
		ptr = kcounts->first;
		while(ptr != NULL) {
			ptr->doubleVal = 0.0;
			ptr = ptr->next;
		}
		
	}
	
	
	// return error
	SEXP rate;
	PROTECT(rate = allocVector(REALSXP, 1));
	REAL(rate)[0] = errcount * 1.0 / len;
	
	dropList(kcounts);
	free(kcounts);
	UNPROTECT(4);
	
	return rate;
	
}		
			
	
	
SEXP mergeClassif(SEXP data, SEXP labels, SEXP KLparam, SEXP rho) {
	// analogous to previous knn, but with vbmerge and building representatives
	
	int i,j, cur, cur2;
	int errcount=0;
	int len = length(labels);
	int levels;
	list *curList;
	SEXP models, plainModel;
	SEXP trueVal, falseVal, ncomp, var;
	PROTECT(trueVal = allocVector(INTSXP, 1));
	INTEGER(trueVal)[0] = 1;
	PROTECT(falseVal = allocVector(INTSXP, 1));
	INTEGER(falseVal)[0] = 0;
	PROTECT(ncomp = allocVector(INTSXP, 1));
	INTEGER(ncomp)[0] = 250;
	//PROTECT(var = allocVector(INTSXP, 1));
	//INTEGER(var)[0] = 128;
	PROTECT(var = allocVector(REALSXP, 1));
	REAL(var)[0] = 0.1;

	double sizes[len];
	
	PROTECT(labels = coerceVector(labels, INTSXP));
	
	if(!isEnvironment(rho)) error("rho should be an environment");
	
	// first regroups item indexes into a metalist. Handle variable structured label vector
	metalist indexes = {NULL, NULL};
	
	for(i=0; i<len; i++) {
		cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		if(cur < 0) {
			addList(&indexes, INTEGER(labels)[i]);
			cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		}
		
		curList = getList(&indexes, cur);
		addItem(curList, i, 0.0);
	}
	
	levels = metaSize(&indexes);
	
	// vector for storing distances
	cell dists[levels];

	// build initial representative vector

	PROTECT(models = allocVector(VECSXP, levels));
	
	Rprintf("initializing controids...\n");
	
	for(i=0; i<levels; i++) {
		PROTECT(plainModel = buildPlainMod(data, listToSXP(getList(&indexes, i)), trueVal));
		SET_VECTOR_ELT(models, i,extractSimpleModel(vbcomp(plainModel, ncomp, var, R_NilValue), falseVal));
		UNPROTECT(1);
	}
	
	// start algorithm
	for(i=0; i<len; i++) {
		
		// first rebuild concerned representative
		cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		curList = getList(&indexes, cur);
		cur2 = indexOfInt(curList, i);
		dropItem(curList, cur2);
		
		Rprintf("merge : computing for item %d, index group %d\n", i, cur);
		
		Rprintf("%d elements merged : ", size(curList));
		for(j=0; j<size(curList); j++) {
			Rprintf("%d ", getIntItem(curList, j));
		}
		Rprintf("\n");
		
		PROTECT(plainModel = buildPlainMod(data, listToSXP(curList), trueVal));
		SET_VECTOR_ELT(models, cur, extractSimpleModel(vbcomp(plainModel, ncomp, var, R_NilValue), falseVal));

		// get size of current centroid for stats
		sizes[i] = length(getListElement(VECTOR_ELT(models, cur), "w"));

		addItem(curList, i, 0.0);
		UNPROTECT(1);
				
		// compute distance of current element wrt representatives
		for(j=0; j<levels; j++) {
			dists[j].ind = getIntMetaItem(&indexes, j);
			dists[j].val = REAL(jsmc(VECTOR_ELT(data, i), VECTOR_ELT(models, j), KLparam))[0];
		}
		
		
		// sort and observe result label : increment error if needed.
		qsort(dists, levels, sizeof(cell), compare);
		
		Rprintf("dists : ");
		for(j=0; j<levels; j++) {
			Rprintf("%f ", dists[j].val);
		}
		Rprintf("\n");

		if(getIntMetaItem(&indexes, cur) != dists[0].ind) {
			errcount++;
		}
		
		Rprintf("seek %d, found %d\n", getIntMetaItem(&indexes, cur), dists[0].ind);
		
		
	}
	
	SEXP res;
	PROTECT(res = allocVector(REALSXP, 2));
	REAL(res)[0] = errcount * 1.0 / len;
	// average size of centroids
	REAL(res)[1] = sum(sizes, len, 1) / len;	
	UNPROTECT(7);
	dropMetalist(&indexes);
	
	return res;
}
		
SEXP constrClassif(SEXP data, SEXP labels, SEXP KLparam, SEXP rho) {
	// analogous to previous knn, but with vbmerge and building representatives
	
	int i,j, cur, cur2;
	int errcount=0;
	int len = length(labels);
	int levels;
	list *curList;
	SEXP models, plainModel;
	SEXP trueVal, falseVal, ncomp, var;
	PROTECT(trueVal = allocVector(INTSXP, 1));
	INTEGER(trueVal)[0] = 1;
	PROTECT(falseVal = allocVector(INTSXP, 1));
	INTEGER(falseVal)[0] = 0;
	PROTECT(ncomp = allocVector(INTSXP, 1));
	INTEGER(ncomp)[0] = 250;
	//PROTECT(var = allocVector(INTSXP, 1));
	//INTEGER(var)[0] = 128;
	PROTECT(var = allocVector(REALSXP, 1));
	REAL(var)[0] = 0.1;

	double sizes[len];

	
	PROTECT(labels = coerceVector(labels, INTSXP));
	
	if(!isEnvironment(rho)) error("rho should be an environment");
	
	// first regroups item indexes into a metalist. Handle variable structured label vector
	metalist indexes = {NULL, NULL};
	
	for(i=0; i<len; i++) {
		cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		if(cur < 0) {
			addList(&indexes, INTEGER(labels)[i]);
			cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		}
		
		curList = getList(&indexes, cur);
		addItem(curList, i, 0.0);
	}
	
	levels = metaSize(&indexes);
	
	// vector for storing distances
	cell dists[levels];

	// build initial representative vector

	PROTECT(models = allocVector(VECSXP, levels));
	
	Rprintf("initializing controids...\n");
	
	for(i=0; i<levels; i++) {
		PROTECT(plainModel = buildPlainMod(data, listToSXP(getList(&indexes, i)), trueVal));
		SET_VECTOR_ELT(models, i,extractSimpleModel(vbconstr(plainModel, ncomp, var, R_NilValue, rho), falseVal));
		UNPROTECT(1);
	}
	
	// start algorithm
	for(i=0; i<len; i++) {
		
		// first rebuild concerned representative
		cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		curList = getList(&indexes, cur);
		cur2 = indexOfInt(curList, i);
		dropItem(curList, cur2);
		
		Rprintf("constr : computing for item %d, index group %d\n", i, cur);
		
		Rprintf("%d elements merged : ", size(curList));
		for(j=0; j<size(curList); j++) {
			Rprintf("%d ", getIntItem(curList, j));
		}
		Rprintf("\n");
		
		PROTECT(plainModel = buildPlainMod(data, listToSXP(curList), trueVal));
		SET_VECTOR_ELT(models, cur, extractSimpleModel(vbconstr(plainModel, ncomp, var, R_NilValue, rho), falseVal));

		// get size of current centroid for stats
		sizes[i] = length(getListElement(VECTOR_ELT(models, cur), "w"));

		addItem(curList, i, 0.0);
		UNPROTECT(1);
				
		// compute distance of current element wrt representatives
		for(j=0; j<levels; j++) {
			dists[j].ind = getIntMetaItem(&indexes, j);
			dists[j].val = REAL(jsmc(VECTOR_ELT(data, i), VECTOR_ELT(models, j), KLparam))[0];
		}
		
		
		// sort and observe result label : increment error if needed.
		qsort(dists, levels, sizeof(cell), compare);
		
		Rprintf("dists : ");
		for(j=0; j<levels; j++) {
			Rprintf("%f ", dists[j].val);
		}
		Rprintf("\n");
		
		if(getIntMetaItem(&indexes, cur) != dists[0].ind) {
			errcount++;
		}
		
		Rprintf("seek %d, found %d\n", getIntMetaItem(&indexes, cur), dists[0].ind);
		
		
	}
	


	SEXP res;
	PROTECT(res = allocVector(REALSXP, 2));
	REAL(res)[0] = errcount * 1.0 / len;
	// average size of centroids
	REAL(res)[1] = sum(sizes, len, 1) / len;

	UNPROTECT(7);
	dropMetalist(&indexes);
	
	return res;
}	


SEXP sampleClassif(SEXP data, SEXP labels, SEXP KLparam, SEXP rho) {
	// analogous to previous knn, but with vbmerge and building representatives
	
	int i,j, cur, cur2;
	int errcount=0;
	int len = length(labels);
	int levels;
	list *curList;
	SEXP models, plainModel;
	SEXP samp;
	SEXP trueVal, falseVal, ncomp, var;
	SEXP nsamp;
	SEXP thres;
	PROTECT(trueVal = allocVector(INTSXP, 1));
	INTEGER(trueVal)[0] = 1;
	PROTECT(falseVal = allocVector(INTSXP, 1));
	INTEGER(falseVal)[0] = 0;
	PROTECT(ncomp = allocVector(INTSXP, 1));
	INTEGER(ncomp)[0] = 250;
	//PROTECT(var = allocVector(INTSXP, 1));
	//INTEGER(var)[0] = 128;
	PROTECT(var = allocVector(REALSXP, 1));
	REAL(var)[0] = 0.1;
	PROTECT(nsamp = allocVector(INTSXP, 1));
	INTEGER(nsamp)[0] = 2000;

	PROTECT(thres = allocVector(REALSXP, 1));
	REAL(thres)[0] = 0.1;

	double sizes[len];

	
	PROTECT(labels = coerceVector(labels, INTSXP));
	
	if(!isEnvironment(rho)) error("rho should be an environment");
	
	// first regroups item indexes into a metalist. Handle variable structured label vector
	metalist indexes = {NULL, NULL};
	
	for(i=0; i<len; i++) {
		cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		if(cur < 0) {
			addList(&indexes, INTEGER(labels)[i]);
			cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		}
		
		curList = getList(&indexes, cur);
		addItem(curList, i, 0.0);
	}
	
	levels = metaSize(&indexes);
	
	// vector for storing distances
	cell dists[levels];

	// build initial representative vector

	PROTECT(models = allocVector(VECSXP, levels));
	
	Rprintf("initializing controids...\n");
	
	for(i=0; i<levels; i++) {
		PROTECT(plainModel = buildPlainMod(data, listToSXP(getList(&indexes, i)), trueVal));
		// Add for sampling
		PROTECT(samp = gmmgen(plainModel, nsamp));
		SET_VECTOR_ELT(models, i,extractSimpleModel(varbayes(VECTOR_ELT(samp, 0), ncomp, thres, R_NilValue), falseVal));
		UNPROTECT(2);
	}
	
	// start algorithm
	for(i=0; i<len; i++) {
		
		// first rebuild concerned representative
		cur = indexOfKey(&indexes, INTEGER(labels)[i]);
		curList = getList(&indexes, cur);
		cur2 = indexOfInt(curList, i);
		dropItem(curList, cur2);
		
		Rprintf("sample : computing for item %d, index group %d\n", i, cur);
		
		Rprintf("%d elements merged : ", size(curList));
		for(j=0; j<size(curList); j++) {
			Rprintf("%d ", getIntItem(curList, j));
		}
		Rprintf("\n");
		

		PROTECT(plainModel = buildPlainMod(data, listToSXP(curList), trueVal));
		PROTECT(samp = gmmgen(plainModel, nsamp));
		SET_VECTOR_ELT(models, cur,extractSimpleModel(varbayes(VECTOR_ELT(samp, 0), ncomp, thres, R_NilValue), falseVal));

		// get size of current centroid for stats
		sizes[i] = length(getListElement(VECTOR_ELT(models, cur), "w"));

		addItem(curList, i, 0.0);
		UNPROTECT(1);
		
		
		// compute distance of current element wrt representatives
		for(j=0; j<levels; j++) {
			dists[j].ind = getIntMetaItem(&indexes, j);
			dists[j].val = REAL(jsmc(VECTOR_ELT(data, i), VECTOR_ELT(models, j), KLparam))[0];
		}
		
		
		// sort and observe result label : increment error if needed.
		qsort(dists, levels, sizeof(cell), compare);
		
		Rprintf("dists : ");
		for(j=0; j<levels; j++) {
			Rprintf("%f ", dists[j].val);
		}
		Rprintf("\n");
		
		if(getIntMetaItem(&indexes, cur) != dists[0].ind) {
			errcount++;
		}
		
		Rprintf("seek %d, found %d\n", getIntMetaItem(&indexes, cur), dists[0].ind);
		
		
	}
	


	SEXP res;
	PROTECT(res = allocVector(REALSXP, 2));
	REAL(res)[0] = errcount * 1.0 / len;
	// average size of centroids
	REAL(res)[1] = sum(sizes, len, 1) / len;

	UNPROTECT(9);
	dropMetalist(&indexes);
	
	return res;
}	

		
		
		
		
		
		
		
