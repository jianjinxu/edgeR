#include "matvec_check.h"

matvec_check::matvec_check(const int nlib, const int nlen, SEXP incoming, const bool transposed, 
		const char* err, const double nullfill) : mycheck(NULL), temp(NULL), isvec(true), istran(transposed), 
		nl(nlib), nt(nlen), tagdex(0), libdex(0) {
	
	std::stringstream failed;
	if (!isNumeric(incoming)) {
		failed << err << " vector or matrix should be double precision";
		throw std::runtime_error(failed.str());
	}
	
	// Checking if it is a vector, matrix or transposed matrix.
	mycheck=REAL(incoming);
	const int curlen=LENGTH(incoming);
	if (curlen==0) {
		temp=new double[nl];
		for (int i=0; i<nl; ++i) { temp[i]=nullfill; }
		mycheck=temp;
	} else if (curlen!=nl) { 
		isvec=false;
 	   	if (LENGTH(incoming)!=nl*nlen) {
			failed << "dimensions of " << err << " matrix are not consistent with number of libraries and tags";
			throw std::runtime_error(failed.str());
		} 
		if (!istran) {
			temp=new double[nl];
			libdex=0;
			for (int i=0; i<nl; ++i, libdex+=nt) { temp[i]=mycheck[libdex]; }
		}
	} else {
		// Otherwise, it's all good; we can use the pointer directly if it's a vector.
		;
	}
	return;
}

void matvec_check::advance() {
	if (!isvec) { 
		if (!istran) { 
			// Copying elements to an array if it is not transposed, so each row can be accessed at a pointer.
			++mycheck;
			if ((++tagdex) >= nt) { return; }
			libdex=0;
			for (int i=0; i<nl; ++i, libdex+=nt) { temp[i]=mycheck[libdex]; }
		} else {
			// Each (original) row is a (transposed) column, so rows can be accessed directly in column-major format.
			mycheck+=nl;
		}
	}
	return;
}

const double* const* matvec_check::access() const { 
	if (isvec || istran) { 
		return &mycheck; 
	} else {
		return &temp;
	}
}

matvec_check::~matvec_check() {
	if (temp!=NULL) { delete [] temp; }
}
