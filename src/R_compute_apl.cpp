#include "glm.h"
#include "matvec_check.h"

SEXP R_compute_apl(SEXP y, SEXP means, SEXP disps, SEXP weights, SEXP adjust, SEXP design) try {
    count_holder counts(y);
    const int num_tags=counts.get_ntags();
    const int num_libs=counts.get_nlibs();
    double* count_ptr=(double*)R_alloc(num_libs, sizeof(double));

    // Setting up the means.
    if (!isReal(means)) {
        throw std::runtime_error("mean matrix must be double-precision");
    }
    if (LENGTH(means)!=LENGTH(y)) {
        throw std::runtime_error("mean and count matrices must be of same size");
    }
    const double* mptr=REAL(means);

    // Setting up the dispersions and weights.
    matvec_check alld(disps, num_tags, num_libs);
    const double* const dptr2=alld.access();
    matvec_check allw(weights, num_tags, num_libs);
    const double* const wptr2=allw.access();

    // Determining whether we want to do the adjustment.
    if (!isLogical(adjust) || LENGTH(adjust)!=1) { 
        throw std::runtime_error("'adjust' must be a logical scalar");
    }
    const bool do_adjust=asLogical(adjust);
    double* W_ptr=(double*)R_alloc(num_libs, sizeof(double));

    // Setting up the design matrix and the CR adjustment object.
    if (!isNumeric(design)) { throw std::runtime_error("design matrix must be double precision"); }
    const int num_coefs=LENGTH(design)/num_libs;
    if (num_coefs*num_libs!=LENGTH(design)) { throw std::runtime_error("dimensions of design matrix not consistent with number of libraries"); }
    adj_coxreid acr(num_libs, num_coefs, REAL(design)); 

    SEXP output=PROTECT(allocVector(REALSXP, num_tags));
    try {
        double* optr=REAL(output);
        double loglik, r, logmur, adj;
        int lib, index;
        for (int tag=0; tag<num_tags; ++tag) {

            /* First computing the log-likelihood. */
            double& sum_loglik=(optr[tag]=0);
            counts.fill_and_next(count_ptr);
            index=0;
            
            for (lib=0; lib<num_libs; ++lib, index+=num_tags) {
                const double& cury=count_ptr[lib];
                const double& curmu=mptr[index];
                if (curmu==0) { 
                    continue; // Should only be possible if count is zero, where the log-likelihood would then be 0.
                }

                if (dptr2[lib] > 0) {
                    // same as loglik <- rowSums(weights*dnbinom(y,size=1/dispersion,mu=mu,log = TRUE))
                    r=1/dptr2[lib];
                    logmur=std::log(curmu+r);
                    loglik = cury*std::log(curmu) - cury*logmur + r*std::log(r) - r*logmur + lgamma(cury+r) - lgamma(cury+1) - lgamma(r); 
                } else {
                    // same as loglik <- rowSums(weights*dpois(y,lambda=mu,log = TRUE))
                    loglik = cury*std::log(curmu) - curmu - lgamma(cury+1);
                }
                sum_loglik += loglik * wptr2[lib]; // with weights.
            }
            
            if (do_adjust) {
                /* Computing 'W', the matrix of negative binomial probabilities. 
                 * The class computes 'XtWX' and performs an LDL decomposition 
                 * to compute the Cox-Reid adjustment factor.
                 */
                for (lib=0; lib<num_libs; ++lib) {
                    const double& curmu=mptr[lib*num_tags];
                    W_ptr[lib] = wptr2[lib] * curmu/(1 + dptr2[lib] * curmu);
                }

                if (num_coefs==1) {
                    adj=0;
                    for (lib=0; lib<num_libs; ++lib) {
                        adj+=W_ptr[lib];
                    }
                    adj=std::log(std::abs(adj))/2;
                } else {
                    std::pair<double, bool> x=acr.compute(W_ptr);
                    if (!x.second) { 
                        std::stringstream errout;
                        errout << "LDL factorization failed for tag " << tag+1;
                        throw std::runtime_error(errout.str());
                    }
                    adj=x.first;
                }
			    sum_loglik-=adj;
			} 

            ++mptr;
            alld.advance();
            allw.advance();
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}

