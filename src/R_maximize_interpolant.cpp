#include "utils.h"
#include "interpolator.h"

SEXP R_maximize_interpolant(SEXP spline_pts, SEXP likelihoods) try {
	if (!isNumeric(spline_pts)) { std::runtime_error("spline points should be a double precision vector"); }
	if (!isNumeric(likelihoods)) { std::runtime_error("likelihoods should be a double precision matrix"); }

    // Loading in the spline x-axis values.
    const int num_pts=LENGTH(spline_pts);
    double* sptr=REAL(spline_pts);
	double* lptr=REAL(likelihoods);
    const int nrows=LENGTH(likelihoods)/num_pts;
    interpolator maxinterpol(num_pts);

   	// Setting up the output object and running through it.
    SEXP output=PROTECT(allocVector(REALSXP, nrows));
    double* out_ptr=REAL(output);
	try { 
		for (int count=0; count<nrows; ++count) {
        	*out_ptr=maxinterpol.find_max(sptr, lptr);
 	   		lptr+=num_pts;	
			++out_ptr;
    	}
	} catch (std::exception &e) {
		UNPROTECT(1);
		throw;
	}
   	UNPROTECT(1);
    return(output);    
} catch (std::exception& e) { return mkString(e.what()); }
