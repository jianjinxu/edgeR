#include "glm.h"
#include "matvec_check.h"

/*** Function to compute the fitted values without a lot of temporary matrices. ***/

SEXP R_get_one_way_fitted (SEXP beta, SEXP offset, SEXP groups) try { 
    SEXP dims=getAttrib(beta, R_DimSymbol);
    if (!isInteger(dims) || LENGTH(dims)!=2) { 
        throw std::runtime_error("matrix dimensions should be an integer vector of length 2");
    }
    int num_tags=INTEGER(dims)[0];
    int num_coefs=INTEGER(dims)[1];

    if (!isReal(beta)) {
        throw std::runtime_error("beta matrix should be double-precision");
    }
    if (LENGTH(beta)!=num_tags*num_coefs) {
        throw std::runtime_error("recorded matrix dimensions are not consistent with its length"); 
    }
    const double* bptr=REAL(beta);
    double* bptr2=(double*)R_alloc(num_coefs, sizeof(double));

    if (!isInteger(groups)) {
        throw std::runtime_error("grouping vector should be integer");
    }
    int num_libs=LENGTH(groups);
    int* gptr=INTEGER(groups);

    matvec_check allo(offset, num_tags, num_libs);
    const double* optr2=allo.access();

    SEXP output=PROTECT(allocMatrix(REALSXP, num_tags, num_libs));
    try {
        double* outptr=REAL(output);
        int lib, coef;
        for (int tag=0; tag<num_tags; ++tag) {
            // Storing to a single vector for faster caching.
            for (coef=0; coef<num_coefs; ++coef) {
                bptr2[coef]=bptr[coef*num_tags];
            }
            ++bptr;

            // Caching is going to be suboptimal, but oh well.
            for (lib=0; lib<num_libs; ++lib) {
                outptr[lib*num_tags]=std::exp(optr2[lib] + bptr2[gptr[lib]]);
            }
            ++outptr;
            allo.advance();
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

