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
    F77_NAME(dsytrf)(&uplo, &ncoefs, working_matrix, &ncoefs, pivots, &temp_work, &lwork, &info);
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
    /* Setting working weight_matrix to 'A=XtWX' with column-major storage.
 	 * This represents the expected Fisher information. The overall strategy is
 	 * to compute the determinant of the information matrix in order to compute
 	 * the adjustment factor for the likelihood (in order to account for uncertainty
 	 * in the nuisance parameters i.e. the fitted values).
 	 */
    for (int row=0; row<ncoefs; ++row) {
        for (int col=0; col<=row; ++col) {
            double& cur_entry=(working_matrix[col*ncoefs+row]=0);
            for (int lib=0; lib<nlibs; ++lib) {
                cur_entry+=design[row*nlibs+lib]*design[col*nlibs+lib]*wptr[lib];
            }
        }
    }   

    /* We now apply the Cholesky decomposition using the appropriate routine from the 
     * LAPACK library. Specifically, we call the routine to do a symmetric indefinite
     * factorisation i.e. A = LDLt. This guarantees factorization for singular matrices 
     * when the actual Cholesky decomposition would fail.
     */
    F77_NAME(dsytrf)(&uplo, &ncoefs, working_matrix, &ncoefs, pivots, work, &lwork, &info);
    if (info<0) { return std::make_pair(0, false); }

    /* For triangular matrices, we need the diagonal to compute the determinant. Fortunately, 
 	 * this remains the diagonal of the output matrix despite the permutations performed by 
 	 * the pivoting. We sum over all log'd diagonal elements to get the log determinant of the 
 	 * information matrix (valid because det(AB)=det(A)*det(B) and we're using the result of 
 	 * the decomposition). We then have to halve the resulting value. 
     */
    double sum_log_diagonals=0;
    for (int i=0; i<ncoefs; ++i) { 
        const double& cur_val=working_matrix[i*ncoefs+i];
        sum_log_diagonals+=(cur_val < low_value ? log_low_value : std::log(cur_val));
    }
	return std::make_pair(sum_log_diagonals*0.5, true);
}


