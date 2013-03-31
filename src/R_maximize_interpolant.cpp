#include "utils.h"
#include "interpolator.h"

extern "C" {

SEXP R_maximize_interpolant(SEXP spline_pts, SEXP likelihoods) try {
	if (!IS_NUMERIC(spline_pts)) { std::runtime_error("spline points should be a double precision vector"); }
	if (!IS_NUMERIC(likelihoods)) { std::runtime_error("likelihoods should be a double precision matrix"); }

    // Loading in the spline x-axis values.
    const int num_pts=LENGTH(spline_pts);
    double* sptr=NUMERIC_POINTER(spline_pts);
	double* lptr=NUMERIC_POINTER(likelihoods);
    const int nrows=LENGTH(likelihoods)/num_pts;
    interpolator maxinterpol(num_pts);

   	// Setting up the output object and running through it.
    SEXP output;
    PROTECT(output=NEW_NUMERIC(nrows));
    double* out_ptr=NUMERIC_POINTER(output);
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

}

