#include "matvec_check.h"

matvec_check::matvec_check(const int nlib, const int nlen, SEXP incoming, const char* err) : nl(nlib) {
	if (!IS_NUMERIC(incoming)) {
		std::stringstream failed;
		failed << err << " matrix should be double precision";
		throw std::runtime_error(failed.str());
	}
	isvec=(LENGTH(incoming)==nl);
	if (!isvec && LENGTH(incoming)!=nl*nlen) {
		std::stringstream failed;
		failed << "dimensions of " << err << " matrix are not consistent with number of libraries and tags";
		throw std::runtime_error(failed.str());
	}
	mycheck=NUMERIC_POINTER(incoming);
	return;
}

void matvec_check::advance() {
	if (!isvec) { mycheck+=nl; }
	return;
}

const double* const matvec_check::access() const { return mycheck; }

