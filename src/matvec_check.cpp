#include "matvec_check.h"

matvec_check::matvec_check(const int nlib, const int nlen, SEXP incoming, const bool transposed, 
		const char* err, const bool nullok) : nl(nlib), nt(nlen), temp(NULL), istran(transposed), 
		isvec(true), libdex(0), tagdex(0) {
	// Checking if NULL (and whether it's allowed). If it is, it becomes a vector of 1's.
	std::stringstream failed;
	if (incoming==R_NilValue) {
		if (!nullok) { 
			failed << err << " vector or matrix cannot be null";
			throw std::runtime_error(failed.str());
		}
		temp=new double[nl];
		for (int i=0; i<nl; ++i) { temp[i]=1; }
		return;
	}
	
	if (!IS_NUMERIC(incoming)) {
		failed << err << " vector or matrix should be double precision";
		throw std::runtime_error(failed.str());
	}
	
	// Checking if it is a vector, matrix or transposed matrix.
	mycheck=NUMERIC_POINTER(incoming);
	if (LENGTH(incoming)!=nl) { 
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
