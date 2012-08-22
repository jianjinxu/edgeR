#include "edgeR.h"

extern "C" {

/*
 * 'w' represents a matrix of negative binomial probabilities (i.e.
 * the prob of success/failure, a function of mean and dispersion)
 * whereas 'x' represents the design matrix. This function calculates
 * the multiplication of the matrices, then performs a Cholesky decomposition
 * to get the lower triangular matrix 'L'.
 */
SEXP cr_adjustment (SEXP w, SEXP x, SEXP nlibs) {
    PROTECT(w=AS_NUMERIC(w));
    PROTECT(x=AS_NUMERIC(x));
    const int num_libs=INTEGER_VALUE(nlibs);
    const long num_tags=LENGTH(w)/((double) num_libs)+0.5;
    const int num_coefs=LENGTH(x)/((double) num_libs)+0.5;
    
    // Setting up a couple of indices for quick access.
    std::deque<double*> w_ptrs(num_libs, NUMERIC_POINTER(w)), 
        x_ptrs(num_coefs, NUMERIC_POINTER(x));
    for (int m=0; m<w_ptrs.size(); ++m) { w_ptrs[m]+=num_tags*m; } 
    for (int n=0; n<x_ptrs.size(); ++n) { x_ptrs[n]+=num_libs*n; } 

    // Setting up variables for Cholesky decomposition. 
    const char uplo='L';                
    int info=0;
    int* pivots=new int[num_coefs];
    double* working_matrix=new double [num_coefs*num_coefs];
    for (int k=0; k<num_coefs*num_coefs; ++k) { working_matrix[k]=0; }

    /* We also want to identify the optimal size of the 'work' array 
     * using the ability of the dystrf function to call ILAENV. We then
     * reallocate the work pointer to this value.
     */
    int lwork=-1;
    double* work=new double[1];
    F77_NAME(dsytrf)(&uplo, &num_coefs, working_matrix, &num_coefs, pivots,
            work, &lwork, &info);
    lwork=(work[0]+0.5);
    delete [] work;
    work=new double[lwork];

    // Running through each tag.
    SEXP output;
    PROTECT(output=allocVector(REALSXP, num_tags));
    double* out_ptr=NUMERIC_POINTER(output);
    for (long tag=0; tag<num_tags; ++tag) {

        // Setting working_matrix to 'A=XWX' with column-major storage.
        for (int row=0; row<num_coefs; ++row) {
            for (int col=0; col<=row; ++col) {
                double& cur_entry=working_matrix[col*num_coefs+row];
                cur_entry=0;
                for (int lib=0; lib<num_libs; ++lib) {
                    cur_entry+=x_ptrs[row][lib]*x_ptrs[col][lib]*w_ptrs[lib][tag];
                }
            }
        }   

        /* We now apply the Cholesky decomposition using the appropriate routine from the 
         * LAPACK library. Specifically, we call the routine to do a symmetric indefinite
         * factorisation i.e. A = LDL*. This guarantees factorization for singular matrices 
         * when the actual Cholesky decomposition would fail.
         */
        F77_NAME(dsytrf)(&uplo, &num_coefs, working_matrix, &num_coefs, pivots, 
                work, &lwork, &info);
        if (info<0) {
            std::stringstream err;
            err << "Factorization failed at line " << tag+1 << ", check count values.";
            error(err.str().c_str());
        }

        /* We now need to extract the diagonal. Fortunately, this remains the diagonal of 
         * the output matrix despite the permutations performed by the pivoting. To get the 
         * Cox-Reid adjustment factor, we sum over all log'd diagonal elements. We then have 
         * to halve the resulting value. 
         */
        double& sum_log_diagonals=(out_ptr[tag]=0);
        for (int i=0; i<num_coefs; ++i) { 
            const double& cur_val=working_matrix[i*num_coefs+i];
            sum_log_diagonals+=(cur_val < low_value ? log_low_value : std::log(cur_val));
        }
        sum_log_diagonals*=0.5;
    } 
                
    UNPROTECT(3);
    delete [] working_matrix;
    delete [] pivots;
    delete [] work;
    return output;
}

}
