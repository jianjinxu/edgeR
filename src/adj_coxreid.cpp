#include "glm.h"
    
const char uplo='L';

adj_coxreid::adj_coxreid (const int& nl, const int& nc, const double* d) : ncoefs(nc), nlibs(nl), info(0), lwork(-1) {
	const int total=ncoefs*ncoefs;
    working_matrix=new double [total];
    for (int k=0; k<total; ++k) { working_matrix[k]=0; }
	pivots=new int [ncoefs];

    /* We also want to identify the optimal size of the 'work' array 
     * using the ability of the dystrf function to call ILAENV. We then
     * reallocate the work pointer to this value.
     */
	double temp_work;
    F77_CALL(dsytrf)(&uplo, &ncoefs, working_matrix, &ncoefs, pivots, &temp_work, &lwork, &info);
	if (info) { 
		delete [] pivots;
		delete [] working_matrix;
		throw std::runtime_error("failed to identify optimal size of workspace through ILAENV"); 
	}
    lwork=int(temp_work+0.5);
    work=new double [lwork];

	// Saving a local copy of the design pointer.
	const int len=nlibs*ncoefs;
	design=new double [len];
	for (int i=0; i<len; ++i) { design[i]=d[i]; }
	return;
}

adj_coxreid::~adj_coxreid () {
	delete [] working_matrix;
	delete [] pivots;
	delete [] work;
	delete [] design;
}

std::pair<double, bool> adj_coxreid::compute(const double* wptr) {
    // Setting working weight_matrix to 'A=Xt %*% diag(W) %*% X' with column-major storage.
    for (int row=0; row<ncoefs; ++row) {
        for (int col=0; col<=row; ++col) {
            double& cur_entry=(working_matrix[col*ncoefs+row]=0);
            for (int lib=0; lib<nlibs; ++lib) {
                cur_entry+=design[row*nlibs+lib]*design[col*nlibs+lib]*wptr[lib];
            }
        }
    }   

    // LDL* decomposition.
    F77_CALL(dsytrf)(&uplo, &ncoefs, working_matrix, &ncoefs, pivots, work, &lwork, &info);
    if (info<0) { return std::make_pair(0, false); }

    // Log-determinant as sum of the log-diagonals, then halving (see below).
    double sum_log_diagonals=0;
    for (int i=0; i<ncoefs; ++i) { 
        const double& cur_val=working_matrix[i*ncoefs+i];
		if (cur_val < low_value || !std::isfinite(cur_val))  { 
			sum_log_diagonals += log_low_value; 
		} else {
			sum_log_diagonals += std::log(cur_val);
		}
    }
	return std::make_pair(sum_log_diagonals*0.5, true);
}

/* EXPLANATION:
   
   XtWX represents the expected Fisher information. The overall strategy is to compute the 
   log-determinant of this matrix, to compute the adjustment factor for the likelihood (in 
   order to account for uncertainty in the nuisance parameters i.e. the fitted values).

   We want to apply the Cholesky decomposition to the XtWX matrix. However, to be safe,
   we call the routine to do a symmetric indefinite factorisation i.e. A = LDLt. This 
   guarantees factorization for singular matrices when the actual Cholesky decomposition 
   would fail because it would start whining about non-positive eigenvectors. 
  
   We then try to compute the determinant of the working_matrix. Here we use two facts:

   - For triangular matrices, the determinant is the product of the diagonals.
   - det(LDL*)=det(L)*det(D)*det(L*)
   - All diagonal elements of 'L' are unity.

   Thus, all we need to do is to we sum over all log'd diagonal elements in 'D', which - 
   happily enough - are stored as the diagonal elements of 'working_matrix'. (And then 
   divide by two, because that's just how the Cox-Reid adjustment works.)

   'info > 0' indicates that one of the diagonals is zero. We handle this by replacing the
   it with an appropriately small non-zero value, if the diagnoal element is zero or NA. This
   is valid because the set of fitted values which are zero will be constant at all dispersions. 
   Thus, any replacement value will eventually cancel out during interpolation to obtain the CRAPLE.

   Note that the LAPACK routine will also do some pivoting, essentially solving PAP* = LDL* for 
   some permutation matrix P. This shouldn't affect anything; the determinant of the permutation 
   is either 1 or -1, but this cancels out, so det(A) = det(PAP*).

   Further note that the routine can theoretically give block diagonals, but this should 
   not occur for positive (semi)definite matrices, which is what XtWX should always be.
*/

