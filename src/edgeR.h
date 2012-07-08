/* Assorted includes and definitions for constants across all the source files
 * for the edgeR shared library.
 */

#ifndef EDGER_H
#define EDGER_H

#include <cmath>
#include <deque>
#include <sstream>
#include <iostream>
#include <algorithm>
extern "C" {
#include "R.h"
#include "Rdefines.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
}

const double low_value=std::pow(10.0, -10.0);
const double log_low_value=std::log(low_value);

#endif
