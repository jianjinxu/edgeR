#include "glm.h"

extern "C" {

SEXP R_levenberg (SEXP nlib, SEXP ntag, SEXP design, SEXP counts, SEXP disp, SEXP offset, SEXP beta, SEXP fitted, SEXP tol, SEXP maxit) try {
	if (!IS_NUMERIC(design)) { throw  std::runtime_error("design matrix should be double precision"); }
	if (!IS_NUMERIC(counts)) { throw  std::runtime_error("count matrix should be double precision"); }
	if (!IS_NUMERIC(disp)) { throw std::runtime_error("dispersion vector should be double precision"); }
	if (!IS_NUMERIC(offset)) { throw std::runtime_error("offset matrix should be double precision"); }
	if (!IS_NUMERIC(beta)) { throw std::runtime_error("matrix of start values for coefficients should be double precision"); }
	if (!IS_NUMERIC(fitted)) { throw std::runtime_error("matrix of starting fitted values should be double precision"); }

    // Getting and checking the dimensions of the arguments.    
    const int num_tags=INTEGER_VALUE(ntag);
    const int num_libs=INTEGER_VALUE(nlib);
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
	} else if (LENGTH(offset)!=clen) {
		throw std::runtime_error("dimensions of offset matrix must match that of the count matrix"); 
	}

    // Initializing pointers to the assorted features.
    double* beta_ptr=NUMERIC_POINTER(beta), *design_ptr=NUMERIC_POINTER(design), *count_ptr=NUMERIC_POINTER(counts), 
		*fitted_ptr=NUMERIC_POINTER(fitted), *offset_ptr=NUMERIC_POINTER(offset), *disp_ptr=NUMERIC_POINTER(disp);

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
			// Copying elements to the new_beta and new_fitted, so output is automatically stored.
			for (int i=0; i<num_libs; ++i) { new_fitted_ptr[i]=fitted_ptr[i]; }
			for (int i=0; i<num_coefs; ++i) { new_beta_ptr[i]=beta_ptr[i]; }
			if (glbg.fit(offset_ptr, count_ptr, *disp_ptr, new_fitted_ptr, new_beta_ptr)) {
				std::stringstream errout;
				errout<< "solution using Cholesky decomposition failed for tag " << tag+1;
				throw std::runtime_error(errout.str());
			} 
			offset_ptr+=num_libs;
			count_ptr+=num_libs;
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

