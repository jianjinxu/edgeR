#include "glm.h"
#include "matvec_check.h"

SEXP R_levenberg (SEXP design, SEXP y, SEXP disp, SEXP offset, SEXP weights,
		SEXP beta, SEXP tol, SEXP maxit) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));
	
    // Getting and checking the dimensions of the arguments.    
	if (!isReal(design)) { throw  std::runtime_error("design matrix should be double precision"); }
    const int dlen=LENGTH(design);
    if (dlen%num_libs!=0) { throw std::runtime_error("size of design matrix is incompatible with number of libraries"); }
    const int num_coefs=dlen/num_libs;
    if (!isReal(beta)) { throw std::runtime_error("starting coefficient values must be positive"); }
    if (LENGTH(beta)!=num_tags*num_coefs) {
        throw std::runtime_error("dimensions of the beta matrix do not match to the number of tags and coefficients");
    } 

    // Initializing pointers to the assorted features.
    const double *design_ptr=REAL(design), *bptr=REAL(beta);
    matvec_check allo(offset, num_tags, num_libs);
    const double* const optr2=allo.access();
    matvec_check allw(weights, num_tags, num_libs);
    const double* const wptr2=allw.access();
    matvec_check alld(disp, num_tags, num_libs);
    const double* const dptr2=alld.access();

    // Initializing output cages.
    SEXP output=PROTECT(allocVector(VECSXP, 5));
   	SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, num_tags, num_coefs)); // beta 
   	SET_VECTOR_ELT(output, 1, allocMatrix(REALSXP, num_tags, num_libs)); // new fitted 
	SET_VECTOR_ELT(output, 2, allocVector(REALSXP, num_tags));
	SET_VECTOR_ELT(output, 3, allocVector(INTSXP, num_tags));
	SET_VECTOR_ELT(output, 4, allocVector(LGLSXP, num_tags));
	double* new_beta_ptr=REAL(VECTOR_ELT(output, 0));
	double* new_fitted_ptr=REAL(VECTOR_ELT(output, 1));
    double* dev_ptr=REAL(VECTOR_ELT(output, 2));
    int* iter_ptr=INTEGER(VECTOR_ELT(output, 3));
    int* fail_ptr=LOGICAL(VECTOR_ELT(output, 4));

    double* tmp_beta_ptr=(double*)R_alloc(num_coefs, sizeof(double));
    double* tmp_fitted_ptr=(double*)R_alloc(num_libs, sizeof(double));

	try {
       	// Running through each tag and fitting the NB GLM.
		glm_levenberg glbg(num_libs, num_coefs, design_ptr, asInteger(maxit), asReal(tol));
        int lib, coef;
    	for (int tag=0; tag<num_tags; ++tag) {
            counts.fill_and_next(count_ptr);
	
			// Copying elements to the tmp_beta as these are modified in-place.
            for (coef=0; coef<num_coefs; ++coef) {
                tmp_beta_ptr[coef]=bptr[coef*num_tags];
            }
            ++bptr;

			if (glbg.fit(optr2, count_ptr, wptr2, dptr2, tmp_fitted_ptr, tmp_beta_ptr)) {
                std::stringstream errout;
				errout<< "solution using Cholesky decomposition failed for tag " << tag+1;
				throw std::runtime_error(errout.str());
			} 
			allo.advance();
			allw.advance();
            alld.advance();
			
            for (lib=0; lib<num_libs; ++lib) {
			    new_fitted_ptr[lib*num_tags]=tmp_fitted_ptr[lib];
            }
            ++new_fitted_ptr;
            for (coef=0; coef<num_coefs; ++coef) {
                new_beta_ptr[coef*num_tags]=tmp_beta_ptr[coef];
            }
            ++new_beta_ptr;

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

