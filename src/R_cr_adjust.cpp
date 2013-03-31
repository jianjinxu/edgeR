#include "glm.h"

extern "C" {

/*
 * 'w' represents a matrix of negative binomial probabilities (i.e.
 * the prob of success/failure, a function of mean and dispersion)
 * whereas 'x' represents the design matrix. This function calculates
 * the multiplication of the matrices, then performs a Cholesky decomposition
 * to get the lower triangular matrix 'L'. The diagonal of 'L' can then
 * be used to compute the Cox-Reid adjustment factor.
 */
SEXP R_cr_adjust (SEXP w, SEXP x, SEXP nlibs) try {
	if (!IS_NUMERIC(w)) { throw std::runtime_error("matrix of likelihoods must be double precision"); }
	if (!IS_NUMERIC(x)) { throw std::runtime_error("design matrix must be double precision"); }

    const int num_libs=INTEGER_VALUE(nlibs);
    const int num_tags=LENGTH(w)/num_libs;
    const int num_coefs=LENGTH(x)/num_libs;
    
    // Setting up a couple of indices for quick access.
    adj_coxreid acr(num_libs, num_coefs, NUMERIC_POINTER(x));
	double* wptr=NUMERIC_POINTER(w);

   	SEXP output=PROTECT(NEW_NUMERIC(num_tags));
    double* out_ptr=NUMERIC_POINTER(output);
    try {
		// Running through each tag.
       	for (long tag=0; tag<num_tags; ++tag) {
			std::pair<double, bool> x=acr.compute(wptr);
			if (!x.second) { 
				std::stringstream errout;
				errout << "LDL factorization failed for tag " << tag+1;
				throw std::runtime_error(errout.str());
			}
			out_ptr[tag]=x.first;
			wptr+=num_libs;
    	} 
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw; 
	}

    UNPROTECT(1);
    return output;
} catch (std::exception& e) { return mkString(e.what()); }

}
