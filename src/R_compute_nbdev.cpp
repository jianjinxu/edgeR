#include "utils.h"
#include "glm.h"
#include "matvec_check.h"

SEXP R_compute_nbdev (SEXP y, SEXP mu, SEXP phi, SEXP weights, SEXP dosum) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));

    // Setting up means.
	if (!isReal(mu)) { throw std::runtime_error("matrix of means should be double-precision"); }
    if (LENGTH(mu)!=num_tags*num_libs) { throw std::runtime_error("length of means is not consistent with length of counts"); }
    const double* mptr=REAL(mu);

    // Setting up dispersions.
    matvec_check alld(phi, num_tags, num_libs);
    const double* const dptr2=alld.access();

    // Seeing if we have to sum things together.
    if (!isLogical(dosum) || LENGTH(dosum)!=1) {
        throw std::runtime_error("summation specification should be a logical vector");
    }
    const bool sumtogether=asLogical(dosum);
    int tag, lib, index;

    if (sumtogether) {
        // Setting up weights.
        matvec_check allw(weights, num_tags, num_libs);
        const double* const wptr2=allw.access();

        SEXP output=PROTECT(allocVector(REALSXP, num_tags));
        try {
            // Running through each row and computing the unit deviance, and then computing the weighted sum.
            double* optr=REAL(output);
            for (tag=0; tag<num_tags; ++tag) {
                counts.fill_and_next(count_ptr);

                double& current_sumdev=(optr[tag]=0);
                index=0;
                for (lib=0; lib<num_libs; ++lib) {
                    current_sumdev += compute_unit_nb_deviance(count_ptr[lib], mptr[index], dptr2[lib]) * wptr2[lib];
                    index+=num_tags;
                }

                ++mptr;
                alld.advance();
                allw.advance();
            }
        } catch (std::exception& e) {
            UNPROTECT(1);
            throw;
        }

        UNPROTECT(1);
        return output;
    } else {
        SEXP output=PROTECT(allocMatrix(REALSXP, num_tags, num_libs));
        try {
            // Computing unit deviance for each observation.
            double* optr=REAL(output);
            for (tag=0; tag<num_tags; ++tag) {
                counts.fill_and_next(count_ptr);

                index=0;
                for (lib=0; lib<num_libs; ++lib) {
                    optr[index] = compute_unit_nb_deviance(count_ptr[lib], mptr[index], dptr2[lib]);
                    index += num_tags;
                }

                ++optr;
                ++mptr;
                alld.advance();
            }
        } catch (std::exception& e) {
            UNPROTECT(1);
            throw;
        }

        UNPROTECT(1);
        return output;
    }
} catch(std::exception& e) {
	return mkString(e.what());
}
