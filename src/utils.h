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

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

const double low_value=std::pow(10.0, -10.0), log_low_value=std::log(low_value);

/* Defining all R-accessible functions. */

extern "C" {

SEXP R_compute_nbdev(SEXP, SEXP, SEXP);

SEXP R_cr_adjust (SEXP, SEXP, SEXP);

SEXP R_exact_test_by_deviance(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_levenberg (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_loess_by_col(SEXP, SEXP, SEXP, SEXP);

SEXP R_maximize_interpolant(SEXP, SEXP);

SEXP R_one_group (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_simple_good_turing (SEXP, SEXP, SEXP);

void processHairpinReads(char**, int*, char**, char**,
		int*, int*, int*, int*, int*, int*,
		int*, int*, int*, int*, char**, int*);
}

#endif
