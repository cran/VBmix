// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <gsl/gsl_math.h>
#include <gsl/gsl_heapsort.h>

#include <math.h>
#include <signal.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "vbconstr.h"
#include "utils.h"
#include "parsers.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>

#define PORT	 1979
#define BUF_SIZE 8192
#define INVALID_SOCKET -1
#define SOCKET_ERROR -1
#define MAXIT 100

// static bounds to index tables: avoid ISO C++ warnings
// ideally, use standard data structures, classes and objects.
#define N_MAX 2000
#define K_MAX 20
#define GROUP_MAX 200


typedef int SOCKET;
typedef struct sockaddr_in SOCKADDR_IN;
typedef struct sockaddr SOCKADDR;
typedef struct in_addr IN_ADDR;

SOCKET sock;

extern "C" {
	SEXP gmmkmsock(SEXP, SEXP, SEXP, SEXP, SEXP);
}

static void interrupt(int val) {
	close(sock);
}


static void init_connection(SEXP host) {
   sock = socket(AF_INET, SOCK_STREAM, 0);
   SOCKADDR_IN sin = { 0 };
   struct hostent *hostinfo;

   if(sock == INVALID_SOCKET) {
      Rprintf("socket()\n");
      exit(errno);
   }

	const char *used_address = CHAR(STRING_ELT(host, 0));
	
   hostinfo = gethostbyname(used_address);
   if (hostinfo == NULL) {
      Rprintf("unknown host %s\n", used_address);
      exit(EXIT_FAILURE);
   }

	Rprintf("connected to %s\n", used_address);

   sin.sin_addr = *(IN_ADDR *) hostinfo->h_addr;
   sin.sin_port = htons(PORT);
   sin.sin_family = AF_INET;

   if(connect(sock,(SOCKADDR *) &sin, sizeof(SOCKADDR)) == SOCKET_ERROR) {
      Rprintf("connect()\n");
      exit(errno);
   }
}


static int compare_doubles (const void * a, const void * b) {
	const double * d_a = (double *)a;
	const double * d_b = (double *)b;
	if (*d_a > *d_b) return 1; 
	else if (*d_a < *d_b) return -1; 
	else 
	return 0; 
}



SEXP gmmkmsock(SEXP models, SEXP names, SEXP ngroups, SEXP rho, SEXP host) {
	// first get info about data input :
	// number of items, dimensionality, number of target groups
	// rho is an environment
	int i, j;
	unsigned found;
	int current;
	
	// we now get the vector of images associated with the data we cluster.
	
	if(!isEnvironment(rho)) error("rho should be an environment");
	
	
	
	const int n = length(models);
	const int k = INTEGER(coerceVector(ngroups, INTSXP))[0];
	const int d = length(VECTOR_ELT(getListElement(VECTOR_ELT(models, 0), "mean"), 0));




	// init connection to server
	// associate specific host
	init_connection(host);
   	char buffer[BUF_SIZE];
   	int nprint;
   	char *ptr;
   	double *doubleptr;
	
	// data structure to monitor set of individuals per centroid (and their distances)
	double dists[N_MAX];
	int inds[N_MAX];
	int tempinds[N_MAX];
	summary summs[K_MAX];
	for(i=0; i<k; i++) {
		summs[i].set = Calloc(n, int);
		summs[i].setdists = Calloc(n, double);
	}
	

	

	// create data structures : 
	// -n length labels vector
	// -k length centroids vector
	
	// we have to manage the case where a centroid might disappear : 
	// prepare to have labels and centroids recast (using reprotect to change protected vector length
	// => additionnally, use an int list to trace effective groups.
	
	// var_ are supposed to remain constant
	SEXP var;
	PROTECT(var = allocVector(REALSXP, 1));
	REAL(var)[0] = 0.1;
	
	SEXP var2;
	PROTECT(var2 = allocVector(INTSXP, 1));
	INTEGER(var2)[0] = 100;
	
	
	SEXP labels;
	PROTECT(labels = allocVector(INTSXP, n));
	// init with it from 1 to n for use with further initialisation
	// we be overriden later
	for(i=0; i<n; i++) {
		INTEGER(labels)[i] = i;
	}
	
	// init groups list
	list* groups = createList();
	for(i=0; i<k; i++) {
		addItem(groups, i, 0.0);
	}
	
	// init centroids size from ngroups
	SEXP centroids;
	PROTECT(centroids = allocVector(VECSXP, k));
	// choose randomly k initial centroids : as executed once we use R function
	SEXP samp;
	SEXP expr, it;

	PROTECT(it = expr = allocVector(LANGSXP, 3));
	
	SETCAR(it, install("sample"));
	it = CDR(it);
	SETCAR(it, labels);
	it = CDR(it);
	SETCAR(it, allocVector(INTSXP, 1));
	INTEGER(CAR(it))[0] = k;
	PROTECT(samp = eval(expr, rho));
	
	// use sample to init centroid structure
	for(i=0; i<k; i++) {
		SET_VECTOR_ELT(centroids, i, VECTOR_ELT(models, INTEGER(samp)[i]));
	}
	
	// immediately unprotect useless structures
	UNPROTECT(2);


	// first protection : memorize protect indexes
	PROTECT_INDEX ipx[3];

	int tmp;
	unsigned cnt = 0;
	unsigned nochange = 0;
	
	// hardcoded MAXIT
	
	while((!nochange) && (cnt < MAXIT)) {
		// first, compute the closest centroid for each item
		// we take care to consider only "active" groups (ie referenced in the list)
		nochange = 1;
		cnt++;
				
		// get changes for output
		list *changes = createList();
		Rprintf("entering div. computing :");
		SEXP curJS;
		PROTECT_WITH_INDEX(curJS=allocVector(REALSXP, 1), &ipx[0]);
		for(i=0; i<n; i++) {
			double cur, curmin = GSL_POSINF;
			int curind=0;
			Rprintf(" %d", i);
			
			// !!! high cost on this loop.
			// weak sample size is used with jsmc => low accuracy.
			// should be optimized (i.e. jsmc that takes set of individuals and set of centroids as param, and factorizes sampling tasks)
			for(j=0; j<size(groups); j++) {
				tmp = getIntItem(groups, j);
				REPROTECT(curJS = jsmc(VECTOR_ELT(models, i), VECTOR_ELT(centroids, tmp), var2), ipx[0]);
				cur = REAL(curJS)[0];
				if(cur<curmin) {
					curind = getIntItem(groups, j);
					curmin = cur;
				}
			}
			
			// keep track of distance of each element wrt its centroid
			dists[i] = curmin;
			
			// control if changed, then affect if necessary
			if(curind != INTEGER(labels)[i]) {
				INTEGER(labels)[i] = curind;
				nochange = 0;
				addItem(changes, i, 0.0);
			}
		}

		UNPROTECT(1);

		Rprintf("\n");
		
		if(!nochange) {
			Rprintf("Following items have changed :");
			for(i=0; i<size(changes); i++) {
				Rprintf(" %d", getIntItem(changes, i));
			}
			Rprintf("\n");
		}
		
		dropList(changes);
		Free(changes);
			
			
		if(!nochange) {
			Rprintf("iteration %d :\n", cnt);
			//printSxpIntVector(labels);
			// intermediary : we analyze current affectations to detect centroids with no associated item.
			for(i=0; i<size(groups); i++) {
				found = 0;
				for(j=0; j<n; j++) {
					if(INTEGER(labels)[j] == getIntItem(groups, i)) {
						found=1;
					}
				}
				if(!found) {
					dropItem(groups, i);
				}
			}
			
			// then, compute new centroids
			SEXP items, curmodels, flag;
			PROTECT_WITH_INDEX(items=allocVector(INTSXP, 1), &ipx[0]);
			PROTECT_WITH_INDEX(flag=allocVector(INTSXP, 1), &ipx[1]);
			PROTECT_WITH_INDEX(curmodels=allocVector(VECSXP, 1), &ipx[2]);
			
			
			for(i=0; i<size(groups); i++) {
				// new list to memorize current set of item indexes
				list *curset = createList();
				for(j=0; j<n; j++) {
					if(INTEGER(labels)[j] == getIntItem(groups, i)) {
						addItem(curset, j, 0.0);
					}
				}
				
				// convert list to a SXP vector
				REPROTECT(items = allocVector(INTSXP, size(curset)), ipx[0]);
				
				
				for(j=0; j<size(curset); j++) {
					INTEGER(items)[j] = getIntItem(curset, j);
				}
				REPROTECT(flag = allocVector(INTSXP, 1), ipx[1]);
				INTEGER(flag)[0] = 1;
				REPROTECT(curmodels = buildPlainMod(models, items, flag), ipx[2]);
				
				// use built model to compute the fusion
				// reuse items to store number of target components : its length is at least 1.
				INTEGER(items)[0] = 150;
				INTEGER(flag)[0] = 0;
				SET_VECTOR_ELT(centroids, getIntItem(groups, i), R_NilValue);
				
				Rprintf("building model %d...\n", getIntItem(groups, i));
				SET_VECTOR_ELT(centroids, getIntItem(groups, i), extractSimpleModel(vbconstr(curmodels, items, var, R_NilValue, rho), flag));
						
				dropList(curset);
				Free(curset);
				//UNPROTECT(3);
			}
			UNPROTECT(3);
			Rprintf("finished computing groups\n");
		}
		
		// we send the information about the clustering structure
		ptr = buffer;
		nprint = sprintf(ptr, "TEMPOSPRING 1.0\n");
		ptr += nprint;
		
		// for each group : build ordered set of indices
		for(i=0; i<size(groups); i++) {
			current = getIntItem(groups, i);
			summs[current].nelts = 0;
			for(j=0; j<n; j++) {
				if(INTEGER(labels)[j] == current) {
					summs[current].set[summs[current].nelts] = j;
					summs[current].setdists[summs[current].nelts] = dists[j];
					summs[current].nelts++;
				}
			}
			// then sort into ascending order
			// inds contains order within setdists
			gsl_heapsort_index((size_t *)inds, summs[current].setdists, summs[current].nelts, sizeof(double), compare_doubles);
			// we keep set apart to perform replacement
			memcpy(tempinds, summs[current].set, n*sizeof(int));
			for(j=0; j<summs[current].nelts; j++) {
				summs[current].set[j] = tempinds[inds[j]];
			}
				
			// write for the current group
			nprint = sprintf(ptr, "%d", current);
			ptr += nprint;
			ptr = writeUrls(summs[current], names, ptr);
			
		}

		// after centroids is built, get distances between them (only for elements present in groups)
		
		Rprintf("computing inter-distances...\n");
		const int groupSize = size(groups);
		double interdists[GROUP_MAX * (GROUP_MAX - 1) / 2];
		doubleptr = interdists;
		for(i=0; i<groupSize-1; i++) {
			for(j=i+1; j<groupSize; j++) {
				*doubleptr = REAL(jsmc(VECTOR_ELT(centroids, getIntItem(groups, i)), VECTOR_ELT(centroids, getIntItem(groups, j)), var2))[0];
				doubleptr++;
			}
		}
		
		*ptr = '\n';
		ptr++;
		
		for(i=0; i<groupSize * (groupSize - 1) / 2; i++) {
			if(interdists[i] == GSL_POSINF) {
				nprint = sprintf(ptr, "9.999e+99");
			} else {
				nprint = sprintf(ptr, "%1.3e", interdists[i]);
			}
			ptr += nprint;
		}
		
		// send string
		*ptr++ = '/';
		*ptr++ = '#';
		*ptr = 0;
		
		send(sock, buffer, strlen(buffer), 0);
		FILE *fp;
		if (fp=fopen("sortie.txt","a"))
		{
			fprintf(fp,"%s\n",buffer);
			fclose(fp);
		}
		
		
		// force garbage collection after each loop
		PROTECT(expr = allocVector(LANGSXP, 1));
	
		SETCAR(expr, install("gc"));
		eval(expr, rho);
		UNPROTECT(1);
			
	}
	
	
	close(sock);
	
	// when set, return structure with "labels" and "centroids"
	// labels is a vector with elements from 1 to k' (for now, we are still from 0 to (k-1) with maybe some missing values)
	// NB : labels is a values vector so no special concern here,
	// but centroids is pointing to input vector... but as long as the original data exists and is not modified integrity is guaranteed.
	// and data is usually immutable.
	
	SEXP targetCent;
	PROTECT(targetCent = allocVector(VECSXP, size(groups)));
	
	// greedy but simple method : get all items in groups, for each reaatribute labels vector and centroid to target.

	for(j=0; j<n; j++) {
		INTEGER(labels)[j] = INTEGER(labels)[j] + 1;
	}
	
	for(i=0; i<n; i++) {
		found = 0;
		j=0;
		while(!found) {
			if(INTEGER(labels)[i] == (getIntItem(groups, j) + 1)) {
				INTEGER(labels)[i] = j+1;
				found = 1;
			} else {
				j++;
			}
		}
	}
				
	
	for(i=0; i<size(groups); i++) {
		SET_VECTOR_ELT(targetCent, i, VECTOR_ELT(centroids, getIntItem(groups, i)));
	}
	
	// create final returned data structure
	SEXP res;
	PROTECT(res = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(res, 0, labels);
	SET_VECTOR_ELT(res, 1, targetCent);
	
	SEXP dims;
	PROTECT(dims = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dims, 0, mkChar("labels"));
	SET_VECTOR_ELT(dims, 1, mkChar("centroids"));

	setAttrib(res, R_NamesSymbol, dims);
	
	// deallocate structures
	
	dropList(groups);
	Free(groups);
	for(i=0; i<k; i++) {
		Free(summs[i].set);
		Free(summs[i].setdists);
	}
	UNPROTECT(6);
	return(res);
}

