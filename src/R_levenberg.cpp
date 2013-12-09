#include "glm.h"
#include "matvec_check.h"

extern "C" {

SEXP R_levenberg (SEXP nlib, SEXP ntag, SEXP design, SEXP counts, SEXP disp, SEXP offset, SEXP weights,
		SEXP beta, SEXP fitted, SEXP tol, SEXP maxit) try {
	if (!IS_NUMERIC(design)) { throw  std::runtime_error("design matrix should be double precision"); }
	if (!IS_NUMERIC(disp)) { throw std::runtime_error("dispersion vector should be double precision"); }
	if (!IS_NUMERIC(beta)) { throw std::runtime_error("matrix of start values for coefficients should be double precision"); }
	if (!IS_NUMERIC(fitted)) { throw std::runtime_error("matrix of starting fitted values should be double precision"); }
    const int num_tags=INTEGER_VALUE(ntag);
    const int num_libs=INTEGER_VALUE(nlib);

	// Checking the count matrix.
    const double *cdptr=NULL;
    const int* ciptr=NULL;
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));
    bool is_integer=IS_INTEGER(counts);
    if (is_integer) {
        ciptr=INTEGER_POINTER(counts);
    } else {
        if (!IS_NUMERIC(counts)) { throw std::runtime_error("count matrix must be integer or double-precision"); }
        cdptr=NUMERIC_POINTER(counts); 
    }
	
    // Getting and checking the dimensions of the arguments.    
    const int dlen=LENGTH(design);
	const int clen=LENGTH(counts);
    if (dlen%num_libs!=0) { throw std::runtime_error("size of design matrix is incompatible with number of libraries"); }
    const int num_coefs=dlen/num_libs;
    if (clen!=num_tags*num_libs) { 
        throw std::runtime_error("dimensions of the count matrix are not as specified");
    } else if (LENGTH(beta)!=num_tags*num_coefs) {
        throw std::runtime_error("dimensions of the beta matrix do not match to the number of tags and coefficients");
    } else if (LENGTH(fitted)!=clen) {
        throw std::runtime_error("dimensions of the fitted matrix do not match those of the count matrix");
    } else if (LENGTH(disp)!=num_tags) { 
		throw std::runtime_error("length of dispersion vector must be equal to the number of tags"); 
	} 

    // Initializing pointers to the assorted features.
    const double* beta_ptr=NUMERIC_POINTER(beta), 
		  *design_ptr=NUMERIC_POINTER(design), 
	  	  *fitted_ptr=NUMERIC_POINTER(fitted), 
		  *disp_ptr=NUMERIC_POINTER(disp);
    matvec_check allo(num_libs, num_tags, offset, true, "offset", false);
    const double* const* optr2=allo.access();
    matvec_check allw(num_libs, num_tags, weights, true, "weight", true);
    const double* const* wptr2=allw.access();

    // Initializing output cages.
    SEXP output=PROTECT(NEW_LIST(5));
   	SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, num_coefs, num_tags)); // beta (transposed)
   	SET_VECTOR_ELT(output, 1, allocMatrix(REALSXP, num_libs, num_tags)); // new fitted (transposed)
	SET_VECTOR_ELT(output, 2, NEW_NUMERIC(num_tags));
	SET_VECTOR_ELT(output, 3, NEW_INTEGER(num_tags));
	SET_VECTOR_ELT(output, 4, NEW_LOGICAL(num_tags));
	double* new_beta_ptr=NUMERIC_POINTER(VECTOR_ELT(output, 0));
	double* new_fitted_ptr=NUMERIC_POINTER(VECTOR_ELT(output, 1));
    double* dev_ptr=NUMERIC_POINTER(VECTOR_ELT(output, 2));
    int* iter_ptr=INTEGER_POINTER(VECTOR_ELT(output, 3));
    int* fail_ptr=LOGICAL_POINTER(VECTOR_ELT(output, 4));
	try {
       	// Running through each tag and fitting the NB GLM.
		glm_levenberg glbg(num_libs, num_coefs, design_ptr, INTEGER_VALUE(maxit), NUMERIC_VALUE(tol));
    	for (int tag=0; tag<num_tags; ++tag) {

			// Copying integer/double counts to a new vector.
            if (is_integer) {
				for (int i=0; i<num_libs; ++i) { count_ptr[i]=double(ciptr[i]); }
				ciptr+=num_libs;
			} else {
				for (int i=0; i<num_libs; ++i) { count_ptr[i]=cdptr[i]; }
				cdptr+=num_libs;
			}
		
			// Copying elements to the new_beta and new_fitted, so output is automatically stored.
			for (int i=0; i<num_libs; ++i) { new_fitted_ptr[i]=fitted_ptr[i]; }
			for (int i=0; i<num_coefs; ++i) { new_beta_ptr[i]=beta_ptr[i]; }
			if (glbg.fit(*optr2, count_ptr, 
#ifdef WEIGHTED
						*wptr2,
#endif
						*disp_ptr, new_fitted_ptr, new_beta_ptr)) {
				std::stringstream errout;
				errout<< "solution using Cholesky decomposition failed for tag " << tag+1;
				throw std::runtime_error(errout.str());
			} 
			allo.advance();
			allw.advance();
			
			++disp_ptr;
			fitted_ptr+=num_libs;
			new_fitted_ptr+=num_libs;
			beta_ptr+=num_coefs;
			new_beta_ptr+=num_coefs;
			
			*(dev_ptr++)=glbg.get_deviance();
			*(iter_ptr++)=glbg.get_iterations();
			*(fail_ptr++)=glbg.is_failure();
    	}
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}

    UNPROTECT(1);
    return output;   
} catch (std::exception& e) { return mkString(e.what()); }

}

