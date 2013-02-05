#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <iostream>
#include <deque>
#include <stdexcept>
#include <sstream>

#include "R.h"
#include "Rdefines.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

const double low_value=std::pow(10.0, -10.0), log_low_value=std::log(low_value);

#endif
