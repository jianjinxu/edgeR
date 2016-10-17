#include "glm.h"
#include "matvec_check.h"

/* Different initialization methods for the Levenberg coefficients */

const char side='L';
const char trans_ormqr='T';
const char uplo='U';
const char trans_trtrs='N';
const char diag='N';
const int unity=1;

struct QRdecomposition {
    QRdecomposition(const double* curX, const int nrows, const int ncoefs) : X(curX), NR(nrows), NC(ncoefs) {
        Xcopy=(double*)R_alloc(NR*NC, sizeof(double));
        tau=(double*)R_alloc(NC, sizeof(double));

        // Also setting up a vector for weights and effects.
        weights=(double*)R_alloc(NR, sizeof(double));
        effects=(double*)R_alloc(NR, sizeof(double));

        // Setting up the workspace for dgeqrf.
        double tmpwork;
        lwork_geqrf=lwork_oqrmqr=-1;
        F77_CALL(dgeqrf)(&NR, &NC, Xcopy, &NR, tau, &tmpwork, &lwork_geqrf, &info);

        // Loading up the optimal WORK.
        lwork_geqrf=tmpwork+0.5;
        work_geqrf=(double*)R_alloc(lwork_geqrf, sizeof(double));

        // Repeating for dormqr
        F77_CALL(dormqr)(&side, &trans_ormqr, &NR, &unity, &NC, Xcopy, &NR, tau, effects, &NR, &tmpwork, &lwork_oqrmqr, &info);
        lwork_oqrmqr=tmpwork+0.5;
        work_oqrmqr=(double*)R_alloc(lwork_oqrmqr, sizeof(double));

        return;
    }

    void store_weights(const double* w) {
        if (w==NULL) {
            std::fill(weights, weights+NR, 1);
        } else {
            for (row=0; row<NR; ++row) {
                weights[row]=std::sqrt(w[row]);
            }
        }
        return;
    }

    void decompose() {
        std::copy(X, X+NR*NC, Xcopy);
        index=0;
        for (coef=0; coef<NC; ++coef) {
            for (row=0; row<NR; ++row) {
                Xcopy[index]*=weights[row];
                ++index;
            }
        }
 
        F77_CALL(dgeqrf)(&NR, &NC, Xcopy, &NR, tau, work_geqrf, &lwork_geqrf, &info);
        if (info) {
            throw std::runtime_error("QR decomposition failed");
        }
       return;
    } 

    void solve(const double * y) {
        for (row=0; row<NR; ++row) {
            effects[row]=y[row]*weights[row];
        }

        F77_CALL(dormqr)(&side, &trans_ormqr, &NR, &unity, &NC, Xcopy, &NR, tau, effects, &NR, work_oqrmqr, &lwork_oqrmqr, &info);
        if (info) {
            throw std::runtime_error("Q**T multiplication failed");
        }

        F77_CALL(dtrtrs)(&uplo, &trans_trtrs, &diag, &NC, &unity, Xcopy, &NR, effects, &NR, &info);
        if (info) {
            throw std::runtime_error("failed to solve the triangular system");
        }

        return;
    }

    const double* X;
    double * Xcopy, * tau, * effects, *weights, *work_geqrf, *work_oqrmqr;
    const int NR, NC;
    int lwork_geqrf, lwork_oqrmqr, info;
    int row, coef, index;
};


SEXP R_get_levenberg_start(SEXP y, SEXP offset, SEXP disp, SEXP weights, SEXP design, SEXP use_null) try {
 	if (!isReal(design)) { throw  std::runtime_error("design matrix should be double precision"); }
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));

    // Getting and checking the dimensions of the arguments.    
    const int dlen=LENGTH(design);
    if (dlen%num_libs!=0) { throw std::runtime_error("size of design matrix is incompatible with number of libraries"); }
    const int num_coefs=dlen/num_libs;

    // Initializing pointers to the assorted features.
    const double *design_ptr=REAL(design);
    matvec_check allo(offset, num_tags, num_libs);
    const double* const optr2=allo.access();
    matvec_check allw(weights, num_tags, num_libs);
    const double* const wptr2=allw.access();
    matvec_check alld(disp, num_tags, num_libs);
    const double* const dptr2=alld.access();

    // Determining what type of algorithm is to be used.
    if (!isLogical(use_null) || LENGTH(use_null)!=1) {
        throw std::runtime_error("'use_null' specification should be a logical scalar");
    }
    const bool null_method=asLogical(use_null);

    SEXP output_beta=PROTECT(allocMatrix(REALSXP, num_tags, num_coefs));
    try {
        double* obptr=REAL(output_beta);
        QRdecomposition QR(design_ptr, num_libs, num_coefs);
        double* tmp_exprs=(double*)R_alloc(num_libs, sizeof(double));
        int lib, coef;

        if (null_method) {
            QR.store_weights(NULL);
            QR.decompose();
            double sum_exprs=0, sum_weight=0, curN, curweight;
            for (int tag=0; tag<num_tags; ++tag) {
                counts.fill_and_next(count_ptr);
               
                // Computing weighted average of the count:library size ratios.
                sum_weight=0; 
                sum_exprs=0;
                for (lib=0; lib<num_libs; ++lib) {
                    curN=std::exp(optr2[lib]);
                    curweight=wptr2[lib]*curN/(1 + dptr2[lib] * curN);
                    sum_exprs += count_ptr[lib] * curweight / curN;
                    sum_weight += curweight;
                }
                std::fill(tmp_exprs, tmp_exprs + num_libs, std::log(sum_exprs/sum_weight));

                // Performing the QR decomposition and taking the solution.
                QR.solve(tmp_exprs);
                for (coef=0; coef<num_coefs; ++coef) {
                    obptr[coef*num_tags]=QR.effects[coef];
                }

                ++obptr;
                allo.advance();
                allw.advance();
                alld.advance();
            }
        } else {
            bool weights_are_the_same=allw.is_row_repeated();
            if (weights_are_the_same) {
                QR.store_weights(wptr2);
                QR.decompose();
            }

            // Finding the delta.
            double delta;
            if (counts.is_data_integer()) {
                const int* ciptr=counts.get_raw_int();
                delta=*std::max_element(ciptr, ciptr+num_libs*num_tags);
            } else {
                const double* cdptr=counts.get_raw_double();
                delta=*std::max_element(cdptr, cdptr+num_libs*num_tags);
            }
            delta=std::min(delta, 1.0/6);

            for (int tag=0; tag<num_tags; ++tag) {
                if (!weights_are_the_same) {
                    QR.store_weights(wptr2);
                    QR.decompose();
                    allw.advance();
                }
                counts.fill_and_next(count_ptr);
              
                // Computing normalized log-expression values.
                for (lib=0; lib<num_libs; ++lib) {
                    tmp_exprs[lib]=std::log(std::max(delta, count_ptr[lib])) - optr2[lib];
                }

                // Performing the QR decomposition and taking the solution.
                QR.solve(tmp_exprs);
                for (coef=0; coef<num_coefs; ++coef) {
                    obptr[coef*num_tags]=QR.effects[coef];
                }

                ++obptr;
                allo.advance();
            }
        }
    } catch (std::exception& e){
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output_beta;
} catch (std::exception& e) {
    return mkString(e.what());
}
