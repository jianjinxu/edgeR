#include "matvec_check.h"

matvec_check::matvec_check(const int nlib, const int nlen, SEXP incoming, const char* err, const bool nullok) : nl(nlib), temp(NULL) {
	std::stringstream failed;
	if (incoming==R_NilValue) {
		if (!nullok) { 
			failed << err << " vector or matrix cannot be null";
			throw std::runtime_error(failed.str());
		}
		temp=new double[nl];
		for (int i=0; i<nl; ++i) { temp[i]=1; }
		mycheck=temp;
		return;
	}
	if (!IS_NUMERIC(incoming)) {
		failed << err << " vector or matrix should be double precision";
		throw std::runtime_error(failed.str());
	}
	isvec=(LENGTH(incoming)==nl);
	if (!isvec && LENGTH(incoming)!=nl*nlen) {
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

matvec_check::~matvec_check() {
	if (temp!=NULL) { delete [] temp; }
}
