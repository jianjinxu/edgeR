#include "matvec_check.h"
#include "utils.h"

bool isNA(int x) {
    return x==NA_INTEGER;
}

bool isNA(double x) {
    return !R_FINITE(x);
}

template <typename T>
SEXP check_counts (const T* ptr, const int len) {
    int allzero=1;
    for (int i=0; i<len; ++i) {
        const T& curval=ptr[i];
        if (isNA(curval)) {
            throw std::runtime_error("missing values not supported");
        }
        if (curval<0) {
            throw std::runtime_error("negative counts not supported");
        }
        if (curval!=0) {
            allzero=0;
        }
    }
    return ScalarLogical(allzero);
}

SEXP R_check_counts(SEXP y) try {
    count_holder counts(y);
    const int total_size=counts.get_ntags()*counts.get_nlibs();
    if (counts.is_data_integer()) {
        return check_counts<int>(counts.get_raw_int(), total_size);
    } else {
        return check_counts<double>(counts.get_raw_double(), total_size);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

SEXP R_check_finite (SEXP values, SEXP name) try {
    if (!isReal(values)) { throw std::runtime_error("should be double-precision"); }
    const int nobs=LENGTH(values);
    const double* vptr=REAL(values);
    for (int o=0; o<nobs; ++o) {
        const double& curval=vptr[o];
        if (!R_FINITE(curval)) {
            throw std::runtime_error("should be finite values");
        }
    }
    return ScalarLogical(1);
} catch (std::exception& e) {
    if (!isString(name) || LENGTH(name)!=1) { throw std::runtime_error("value-type name should be a string"); }
    std::stringstream final;
    final << CHAR(STRING_ELT(name, 0)) << " " << e.what();
    return mkString(final.str().c_str());
}

SEXP R_check_positive (SEXP values, SEXP name) try {
    if (!isReal(values)) { throw std::runtime_error("should be double-precision"); }
    const int nobs=LENGTH(values);
    const double* vptr=REAL(values);
    for (int o=0; o<nobs; ++o) {
        const double& curval=vptr[o];
        if (!R_FINITE(curval) || curval <= 0) {
            throw std::runtime_error("should be finite positive values");
        }
    }
    return ScalarLogical(1);
} catch (std::exception& e) {
     if (!isString(name) || LENGTH(name)!=1) { throw std::runtime_error("value-type name should be a string"); }
    std::stringstream final;
    final << CHAR(STRING_ELT(name, 0)) << " " << e.what();
    return mkString(final.str().c_str());
}

SEXP R_check_nonnegative (SEXP values, SEXP name) try {
    if (!isReal(values)) { throw std::runtime_error("should be double-precision"); }
    const int nobs=LENGTH(values);
    const double* vptr=REAL(values);
    for (int o=0; o<nobs; ++o) {
        const double& curval=vptr[o];
        if (!R_FINITE(curval) || curval < 0) {
            throw std::runtime_error("should be finite non-negative values");
        }
    }
    return ScalarLogical(1);
} catch (std::exception& e) {
     if (!isString(name) || LENGTH(name)!=1) { throw std::runtime_error("value-type name should be a string"); }
    std::stringstream final;
    final << CHAR(STRING_ELT(name, 0)) << " " << e.what();
    return mkString(final.str().c_str());
}

template <typename T>
SEXP check_zero_fitted(const T* yptr, const int num_tags, const int num_libs, SEXP fitted, SEXP tolerance) {
    const int total_len=num_tags*num_libs;

    if (!isReal(fitted)) { throw std::runtime_error("fitted values must be double-precision"); }
    if (LENGTH(fitted)!=num_tags*num_libs) { throw std::runtime_error("dimensions of fitted and count matrices are not equal"); }
    const double* fptr=REAL(fitted);

    if (!isReal(tolerance) || LENGTH(tolerance)!=1) { throw std::runtime_error("tolerance must be a double-precision scalar"); }
    const double tol=asReal(tolerance);

    SEXP output=PROTECT(allocMatrix(LGLSXP, num_tags, num_libs));
    try {
        int* optr=LOGICAL(output);
        for (int i=0; i<total_len; ++i) {
            (*optr)=((*fptr < tol) & (*yptr < tol));
            ++optr;
            ++fptr;
            ++yptr;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP R_check_zero_fitted(SEXP y, SEXP fitted, SEXP tolerance) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();

    if (counts.is_data_integer()){ 
        return check_zero_fitted<int>(counts.get_raw_int(), num_tags, num_libs, fitted, tolerance);
    } else {
        return check_zero_fitted<double>(counts.get_raw_double(), num_tags, num_libs, fitted, tolerance);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}
