// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include "parsers.h"
#define BUF_SIZE 128

void transformPath(char *chaine) {
	int i;
	char *cur = chaine;
	for(i=0; i<strlen(chaine); i++) {
		if(*cur == '/') {
			memmove(cur+1, cur, strlen(cur)+1);
			*cur = '\\';
			cur++;
			*cur = '\\';
		}
		cur++;
	}
}
		

char *writeUrls(summary summ, SEXP names, char *buffer) {
	// writes urls into buffer, we return the pointer to continue writing in buffer
	
	// first decide how many urls will be written : depends on summ->nelts value
	int leftw;
	int rightw;
	int nprint;
	int i;
	char *ptr = buffer;
	char current[BUF_SIZE];
	int currentIndex;
	
	int n = summ.nelts;
	int temp;
	temp = n % 2;
	
	if(n > 7) {
		leftw = 4;
		rightw=4;
	} else if(n == 1) {
		leftw=1;
		rightw=0;
	} else if(temp == 0) {
		leftw = n / 2;
		rightw = n / 2;
	} else {
		leftw = (n-1)/2;
		rightw=(n-1)/2;
	}
	
	// then right appropriate output
	for(i=0; i<leftw; i++) {
		if(i == 0) {
			nprint = sprintf(ptr, "/");
			ptr += nprint;
		} else {
			nprint = sprintf(ptr, "#");
			ptr += nprint;
		}
		currentIndex = summ.set[i];
		strcpy(current, CHAR(STRING_ELT(names, currentIndex)));
		transformPath(current);
		nprint = sprintf(ptr, "%s", current);
		ptr += nprint;
	}
	
	*ptr = '/';
	ptr++;
	
	for(i=0; i<rightw; i++) {
		if(i != 0) {
			nprint = sprintf(ptr, "#");
			ptr += nprint;
		}
		currentIndex = summ.set[n-1-i];
		strcpy(current, CHAR(STRING_ELT(names, currentIndex)));
		transformPath(current);
		nprint = sprintf(ptr, "%s", current);
		ptr += nprint;
	}
	
	*ptr = '/';
	ptr++;
	
	return(ptr);		
		
}






