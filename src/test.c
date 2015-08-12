
#ifndef TEST_C
#define TEST_C

#include "test.h"

SEXP test(SEXP arg) {
	int c_arg = INTEGER(arg)[0];
	if(c_arg == TRUE) {
		Rprintf("TRUE detected");
	} else {
		Rprintf("FALSE detected");
	}

	return(R_NilValue);

}

#endif
