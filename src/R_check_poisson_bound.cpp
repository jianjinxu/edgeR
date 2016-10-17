#include "matvec_check.h"

/* Checks whether the variance is below the Poisson bound. */

SEXP R_check_poisson_bound (SEXP fitted, SEXP disp, SEXP s2) try {
    if (!isReal(fitted)) { throw std::runtime_error("matrix of fitted values should be double-precision"); }
    const double* fptr=REAL(fitted);

    SEXP dims=getAttrib(fitted, R_DimSymbol);
    if (!isInteger(dims) || LENGTH(dims)!=2) { 
        throw std::runtime_error("matrix dimensions should be an integer vector of length 2");
    }
    const int num_tags=INTEGER(dims)[0], num_libs=INTEGER(dims)[1];
    if (LENGTH(fitted)!=num_libs*num_tags) {
        throw std::runtime_error("recorded matrix dimensions are not consistent with its length"); 
    }

    matvec_check alld(disp, num_tags, num_libs);
    const double* const dptr2=alld.access();
    matvec_check alls(s2, num_tags, num_libs);
    const double* const sptr2=alls.access();

    SEXP output=PROTECT(allocVector(LGLSXP, num_tags));
    try {
        int* optr=LOGICAL(output);
        std::fill(optr, optr+num_tags, 0);
        
        int lib;
        for (int tag=0; tag<num_tags; ++tag) {
            for (lib=0; lib<num_libs; ++lib) {
                if ((fptr[lib*num_tags] * dptr2[lib] + 1) * sptr2[lib] < 1) {
                    optr[tag]=1;
                    break;
                }
            }

            ++fptr;
            alld.advance();
            alls.advance();
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e){ 
    return mkString(e.what());
}
