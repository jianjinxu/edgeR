#include "glm.h"
#include "add_prior.h"
#include "matvec_check.h"

SEXP R_ave_log_cpm(SEXP y, SEXP offset, SEXP prior, SEXP disp, SEXP weights, SEXP max_iterations, SEXP tolerance) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));

    add_prior AP(num_tags, num_libs, prior, offset, true, true);
    const double* const out_prior=AP.get_priors();
    const double* const out_off=AP.get_offsets();
    const bool priors_are_the_same=AP.same_across_rows();

    // Other CompressedMatrix stuff.
    matvec_check alld(disp, num_tags, num_libs);
    const double* const dptr2=alld.access();
    matvec_check allw(weights, num_tags, num_libs);
    const double* const wptr2=allw.access();
    
    // GLM fitting specifications
    const int maxit=asInteger(max_iterations);
    const double tol=asReal(tolerance);

    SEXP output=PROTECT(allocVector(REALSXP, num_tags));
    try {
        double* optr=REAL(output);
        if (priors_are_the_same) {
            AP.fill_and_next();
        }
        
        int lib;
        for (int tag=0; tag<num_tags; ++tag) {
            counts.fill_and_next(count_ptr);
            if (!priors_are_the_same) {
                AP.fill_and_next();
            }

            // Adding the current set of priors.
            for (lib=0; lib<num_libs; ++lib) {
                count_ptr[lib]+=out_prior[lib];
            }

            // Fitting a one-way layout.
            std::pair<double,bool> fit=glm_one_group(num_libs, maxit, tol, out_off, wptr2, count_ptr, dptr2, NA_REAL);
            optr[tag]=(fit.first + LNmillion)/LNtwo;

            allw.advance();
            alld.advance();
        }

    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}
