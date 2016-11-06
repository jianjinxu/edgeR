#include "matvec_check.h"
#include "utils.h"
#include "add_prior.h"

/**** Calculating the CPMs in cpm.default with log=TRUE, but more memory-efficiently. *******/

SEXP R_calculate_cpm_log (SEXP y, SEXP libsize, SEXP prior) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));

    add_prior AP(num_tags, num_libs, prior, libsize, false, true);
    const double* const out_prior=AP.get_priors();
    const double* const out_off=AP.get_offsets();
    const bool priors_are_the_same=AP.same_across_rows();

    SEXP output=PROTECT(allocMatrix(REALSXP, num_tags, num_libs));
    try {
        double* outptr=REAL(output);
        if (priors_are_the_same) { // Just doing this once to save time, if they're all the same.
            AP.fill_and_next();
        }

        // Adding a library size-adjusted prior to each count.
        int lib;
        for (int tag=0; tag<num_tags; ++tag) { 
            counts.fill_and_next(count_ptr);

            if (!priors_are_the_same) { // Repeating with the next set of priors/offsets.
                AP.fill_and_next();
            }

            for (lib=0; lib<num_libs; ++lib) {
                double& curval=outptr[lib*num_tags];
                curval=count_ptr[lib] + out_prior[lib];
                curval=std::log(curval)-out_off[lib]+LNmillion;
                curval/=LNtwo;
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

/**** Calculating the CPMs in cpm.default with log=FALSE, but more memory-efficiently. *******/

SEXP R_calculate_cpm_raw (SEXP y, SEXP libsize) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));
    matvec_check allL(libsize, num_tags, num_libs);
    const double* const lptr2=allL.access();

    SEXP output=PROTECT(allocMatrix(REALSXP, num_tags, num_libs));
    try {
        double* outptr=REAL(output);
        int lib;
        for (int tag=0; tag<num_tags; ++tag) { 
            counts.fill_and_next(count_ptr);

            for (lib=0; lib<num_libs; ++lib) {
                const double& curlib=lptr2[lib];
                outptr[lib*num_tags]=(curlib==0 ? R_NaN : count_ptr[lib]/(curlib/one_million));
            }
            ++outptr;
            allL.advance();
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



