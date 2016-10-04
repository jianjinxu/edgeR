#include "matvec_check.h"
#include "utils.h"
#include "add_prior.h"

/**** Adding a prior count to each observation. *******/

SEXP R_add_prior_count(SEXP y, SEXP offset, SEXP prior) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));

    add_prior AP(num_tags, num_libs, prior, offset, true, true);
    const double* const out_prior=AP.get_priors();
    const double* const out_off=AP.get_offsets();
    const bool priors_are_the_same=AP.same_across_rows();

    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, num_tags, num_libs));
        double* outptr=REAL(VECTOR_ELT(output, 0));

        double* libptr;
        if (priors_are_the_same) {
            // Just doing this once to save time, if they are all the same.
            AP.fill_and_next();
            SET_VECTOR_ELT(output, 1, allocVector(REALSXP, num_libs));
            libptr=REAL(VECTOR_ELT(output, 1));
            std::copy(out_off, out_off+num_libs, libptr);
        } else {
            SET_VECTOR_ELT(output, 1, allocMatrix(REALSXP, num_tags, num_libs));
            libptr=REAL(VECTOR_ELT(output, 1));
        }

        // Adding a library size-adjusted prior to each count.
        int lib;
        for (int tag=0; tag<num_tags; ++tag) { 
            counts.fill_and_next(count_ptr);

            if (!priors_are_the_same) {
                // Repeating with the next set of priors/offsets, and storing the new offsets.
                AP.fill_and_next();
                for (lib=0; lib<num_libs; ++lib) {
                    libptr[lib*num_tags]=out_off[lib];
                }
                ++libptr;
            }

            for (lib=0; lib<num_libs; ++lib) {
                outptr[lib*num_tags]=count_ptr[lib] + out_prior[lib];
            }
            ++outptr;
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

