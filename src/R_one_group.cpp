#include "glm.h"
#include "matvec_check.h"

template <typename T>
bool is_array_equal_to(const T* x, const int n, const bool rep, const T& v) {
    if (rep) {
        return (n>0 && x[0]==v);
    } else {
        for (int i=0; i<n; ++i) {
            if (x[i]!=v) { return false; }
        }
        return true;
    }
}

SEXP R_one_group (SEXP y, SEXP disp, SEXP offsets, SEXP weights, SEXP max_iterations, SEXP tolerance, SEXP beta) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
  	double* yptr=(double*)R_alloc(num_libs, sizeof(double));
 
    // Setting up assorted input matrices.
    matvec_check allo(offsets, num_tags, num_libs);
	const double* const optr2=allo.access();
	matvec_check allw(weights, num_tags, num_libs);
	const double* const wptr2=allw.access();
    matvec_check alld(disp, num_tags, num_libs);
	const double* const dptr2=alld.access();
    matvec_check allb(beta, num_tags, 1); // only one coefficient.
    const double* const bptr2=allb.access();
    
    // GLM iterations.
	const int maxit=asInteger(max_iterations);
	const double tol=asReal(tolerance);
   
    // Setting up beta for output. 
	SEXP output=PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(output, 0, allocVector(REALSXP, num_tags));
	SET_VECTOR_ELT(output, 1, allocVector(LGLSXP, num_tags));
    double* bptr=REAL(VECTOR_ELT(output, 0));
	int* cptr=LOGICAL(VECTOR_ELT(output, 1));
	try {

        // Preparing for possible Poisson sums.
        bool disp_is_zero, weight_is_one;
        double sum_counts, sum_lib=0;
        int lib;
        if (allo.is_row_repeated()) { 
            for (lib=0; lib<num_libs; ++lib) {
                sum_lib+=std::exp(optr2[lib]);
            }
        }
        if (alld.is_row_repeated()) {
            disp_is_zero=is_array_equal_to<double>(dptr2, num_libs, alld.is_col_repeated(), 0);
        }
        if (allw.is_row_repeated()) {
            weight_is_one=is_array_equal_to<double>(wptr2, num_libs, allw.is_col_repeated(), 1);
        }
        
    	// Iterating through tags and fitting.
    	for (int tag=0; tag<num_tags; ++tag) {
            counts.fill_and_next(yptr);
    
            // Checking for the Poisson special case with all-unity weights and all-zero dispersions.
            if (!alld.is_row_repeated()) {
                disp_is_zero=is_array_equal_to<double>(dptr2, num_libs, alld.is_col_repeated(), 0);
            }
            if (!allw.is_row_repeated()) {
                weight_is_one=is_array_equal_to<double>(wptr2, num_libs, allw.is_col_repeated(), 1);
            }

            if (disp_is_zero && weight_is_one) {
                if (!allo.is_row_repeated()) {
                    // Only recalculate sum of library sizes if it has changed.
                    sum_lib=0;
                    for (lib=0; lib<num_libs; ++lib) { sum_lib+=std::exp(optr2[lib]); }
                }

                sum_counts=0;
                for (lib=0; lib<num_libs; ++lib) { sum_counts+=yptr[lib]; }
                if (sum_counts==0) {
                    bptr[tag]=R_NegInf;
                } else {
                    bptr[tag]=std::log(sum_counts/sum_lib);
                }
                cptr[tag]=true;
            } else {
                // Otherwise going through NR iterations.
                std::pair<double, bool> out=glm_one_group(num_libs, maxit, tol, optr2, wptr2, yptr, dptr2, *bptr2);
                bptr[tag]=out.first;
                cptr[tag]=out.second;
            }

			allo.advance();
			allw.advance();
            alld.advance();
            allb.advance();
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
