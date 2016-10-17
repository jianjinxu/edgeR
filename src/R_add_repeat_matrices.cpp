#include "matvec_check.h"
#include "utils.h"

SEXP R_add_repeat_matrices(SEXP x, SEXP y, SEXP nr, SEXP nc) try {
    if (!isInteger(nr) || LENGTH(nr)!=1) { throw std::runtime_error("number of rows should be an integer scalar"); }
    const int nrows=asInteger(nr);
    if (!isInteger(nc) || LENGTH(nc)!=1) { throw std::runtime_error("number of columns should be an integer scalar"); }
    const int ncols=asInteger(nc);

    matvec_check allx(x, nrows, ncols);
    const double* const xptr2=allx.access();
    matvec_check ally(y, nrows, ncols);
    const double* const yptr2=ally.access();

    SEXP output=PROTECT(allocVector(VECSXP, 3));
    try {
        SET_VECTOR_ELT(output, 0, allocMatrix(REALSXP, nrows, ncols));
        double* optr=REAL(VECTOR_ELT(output, 0));
        int tag, lib;

        for (tag=0; tag<nrows; ++tag) {
            for (lib=0; lib<ncols; ++lib) {
                optr[lib*nrows] = xptr2[lib] + yptr2[lib];
            }
            ++optr;
            allx.advance();
            ally.advance();
        }

        SET_VECTOR_ELT(output, 1, ScalarLogical(allx.is_row_repeated() & ally.is_row_repeated()));
        SET_VECTOR_ELT(output, 2, ScalarLogical(allx.is_col_repeated() & ally.is_col_repeated()));
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}
