#ifndef UTILS_H
#define UTILS_H
//#define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

#include <cmath>
#include <deque>
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

/* Defining all R-accessible functions. */

extern "C" {

SEXP R_compute_nbdev(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_compute_apl (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_exact_test_by_deviance(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_levenberg (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_get_levenberg_start (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_loess_by_col(SEXP, SEXP, SEXP, SEXP);

SEXP R_maximize_interpolant(SEXP, SEXP);

SEXP R_one_group (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_get_one_way_fitted (SEXP, SEXP, SEXP);

SEXP R_simple_good_turing (SEXP, SEXP, SEXP);

SEXP R_add_prior_count (SEXP, SEXP, SEXP);

SEXP R_calculate_cpm_log (SEXP, SEXP, SEXP);

SEXP R_calculate_cpm_raw (SEXP, SEXP);

SEXP R_ave_log_cpm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_check_counts (SEXP);

SEXP R_check_finite (SEXP, SEXP);

SEXP R_check_positive (SEXP, SEXP);

SEXP R_check_nonnegative (SEXP, SEXP);

SEXP R_check_zero_fitted (SEXP, SEXP, SEXP);

SEXP R_check_poisson_bound (SEXP, SEXP, SEXP);

SEXP R_add_repeat_matrices(SEXP, SEXP, SEXP, SEXP);

void processHairpinReads(int *, int *, char**, char**, int*,
		char**, char**, int*, int*, int*, int*, int*, int*,
		int*, int*, int*, int*, int*, int*, int *,
		int *, char**, int*);

}

/* Other utility functions and values */

const double low_value=std::pow(10.0, -10.0), log_low_value=std::log(low_value);

const double LNtwo=std::log(2), one_million=1000000, LNmillion=std::log(one_million);

#endif
