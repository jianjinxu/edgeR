#include "glm.h"

extern "C" {

SEXP R_one_group (SEXP nt, SEXP nl, SEXP y, SEXP disp, SEXP offsets, SEXP max_iterations, SEXP tolerance) try {
	const int num_tags=INTEGER_VALUE(nt);
	const int num_libs=INTEGER_VALUE(nl);
	if (num_tags*num_libs != LENGTH(y) ) { throw std::runtime_error("dimensions of the count table are not as specified"); }  // Checking that it is an exact division.
  
	const int maxit=INTEGER_VALUE(max_iterations);
	const double tol=NUMERIC_VALUE(tolerance);
	if (!IS_NUMERIC(disp)) { throw std::runtime_error("dispersion vector must be double precision"); }
	if (!IS_NUMERIC(offsets))  { throw std::runtime_error("offset matrix/vector must be double precision"); }

	if (LENGTH(disp)!=num_tags) { throw std::runtime_error("length of dispersion vector must be 1 or equal to the number of tags"); }
   	if (LENGTH(offsets)!=num_tags*num_libs) { throw std::runtime_error("dimensions of offset matrix must match that of the count matrix"); }
    
    // Setting up some iterators. We provide some flexibility to detecting numeric-ness.
	bool is_integer=IS_INTEGER(y);
	double *ydptr=0;
	int* yiptr=0;
	if (is_integer) { 
		yiptr=INTEGER_POINTER(y); 
		ydptr=(double*) R_alloc(num_libs, sizeof(double));
	} else { 
		ydptr=NUMERIC_POINTER(y); 
	}
	double* optr=NUMERIC_POINTER(offsets);
	double* dptr=NUMERIC_POINTER(disp);

    // Setting up beta for output. 
	SEXP output=PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(output, 0, NEW_NUMERIC(num_tags));
	SET_VECTOR_ELT(output, 1, NEW_LOGICAL(num_tags));
    double* bptr=NUMERIC_POINTER(VECTOR_ELT(output, 0));
	int* cptr=LOGICAL_POINTER(VECTOR_ELT(output, 1));
	try {
        
    	// Iterating through tags and fitting.
    	for (int tag=0; tag<num_tags; ++tag) {
			if (is_integer) { 
				for (int i=0; i<num_libs; ++i) { ydptr[i]=yiptr[i]; }	
				yiptr+=num_libs;
			}
			std::pair<double, bool> out=glm_one_group(num_libs, maxit, tol, optr, ydptr, *dptr);
			(*bptr)=out.first;
			(*cptr)=out.second;
			if (!is_integer) { ydptr+=num_libs; }
			optr+=num_libs;
			++bptr;
			++cptr;
			++dptr;
    	}
	} catch (std::exception& e) { 
		UNPROTECT(1);
		throw; 
	}

	// Returning everything as a list.
    UNPROTECT(1); 
    return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

}
