#include "utils.h"
#include "glm.h"

SEXP R_compute_nbdev (SEXP y, SEXP mu, SEXP phi) try {
	if (!isNumeric(phi)) { throw std::runtime_error("dispersion vector should be double-precision"); }
	const int ntags=LENGTH(phi);
	if (!isNumeric(y)) { throw std::runtime_error("count matrix should be double-precision"); }
	if (!isNumeric(mu)) { throw std::runtime_error("matrix of means should be double-precision"); }
	const int nlib=LENGTH(mu)/ntags;
	if (nlib*ntags !=LENGTH(mu)) { throw std::runtime_error("mean matrix has inconsistent dimensions"); }
	if (LENGTH(mu)!=LENGTH(y)) { throw std::runtime_error("count and mean matrices should have same dimensions"); }

	const double* yptr=REAL(y);
	const double* mptr=REAL(mu);
	const double* dptr=REAL(phi);

	// Running through each row and computing the unit deviance, and then that sum.
	SEXP output=PROTECT(allocMatrix(REALSXP, ntags, nlib));
	try {
		double* optr=REAL(output);
		int counter;
		for (int i=0; i<ntags; ++i) {
			counter=0;
			for (int j=0; j<nlib; ++j, counter+=ntags) {
				optr[counter]=compute_unit_nb_deviance(yptr[counter], mptr[counter], dptr[i]);
			}
			++optr;
			++yptr;
			++mptr;
		}
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}

	UNPROTECT(1);
	return output;
} catch(std::exception& e) {
	return mkString(e.what());
}
