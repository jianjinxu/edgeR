#include "glm.h"
#include "matvec_check.h"

SEXP R_one_group (SEXP nt, SEXP nl, SEXP y, SEXP disp, SEXP offsets, SEXP weights, SEXP max_iterations, SEXP tolerance) try {
	const int num_tags=asInteger(nt);
	const int num_libs=asInteger(nl);
	if (num_tags*num_libs != LENGTH(y) ) { throw std::runtime_error("dimensions of the count table are not as specified"); }  // Checking that it is an exact division.
  
	const int maxit=asInteger(max_iterations);
	const double tol=asReal(tolerance);
	if (!isNumeric(disp)) { throw std::runtime_error("dispersion vector must be double precision"); }
	if (LENGTH(disp)!=num_tags) { throw std::runtime_error("length of dispersion vector must be 1 or equal to the number of tags"); }
    
    // Setting up some iterators. We provide some flexibility to detecting numeric-ness.
	double *ydptr=NULL;
	int* yiptr=NULL;
	double* yptr=(double*)R_alloc(num_libs, sizeof(double));
	bool is_integer=isInteger(y);
	if (is_integer) { 
		yiptr=INTEGER(y); 
	} else { 
		if (!isNumeric(y)) { throw std::runtime_error("count matrix must be integer or double-precision"); }
		ydptr=REAL(y); 
	}
    matvec_check allo(num_libs, num_tags, offsets, false, "offset", false);
	const double* const* optr2=allo.access();
	matvec_check allw(num_libs, num_tags, weights, false, "weight", true);
	const double* const* wptr2=allw.access();
	double* dptr=REAL(disp);

    // Setting up beta for output. 
	SEXP output=PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(output, 0, allocVector(REALSXP, num_tags));
	SET_VECTOR_ELT(output, 1, allocVector(LGLSXP, num_tags));
    double* bptr=REAL(VECTOR_ELT(output, 0));
	int* cptr=LOGICAL(VECTOR_ELT(output, 1));
	try {
        
    	// Iterating through tags and fitting.
    	int counter=0;
    	for (int tag=0; tag<num_tags; ++tag) {
			counter=0;
			if (is_integer) { 
				for (int i=0; i<num_libs; ++i, counter+=num_tags) { yptr[i]=double(yiptr[counter]); }	
				++yiptr;
			} else {
				for (int i=0; i<num_libs; ++i, counter+=num_tags) { yptr[i]=ydptr[counter]; }
				++ydptr;
			}
			std::pair<double, bool> out=glm_one_group(num_libs, maxit, tol, *optr2,
#ifdef WEIGHTED					
					*wptr2,
#endif					
					yptr, *dptr);

			(*bptr)=out.first;
			(*cptr)=out.second;
			++bptr;
			++cptr;
			++dptr;
			allo.advance();
			allw.advance();
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
