#include "edgeR.h"

/* Function to calculate the deviance. Note the protection for very large mu*phi (where we
 * use a gamma instead) or very small mu*phi (where we use the Poisson instead). This 
 * approximation protects against numerical instability introduced by subtrackting
 * a very large log value in (log cur_mu) with another very large logarithm (log cur_mu+1/phi).
 * We need to consider the 'phi' as the approximation is only good when the product is
 * very big or very small.
 */
    
const double one_million=std::pow(10, 6.0), one_millionth=std::pow(10, -6.0),
     mildly_low_value=std::pow(10, -8.0), supremely_low_value=std::pow(10, -13.0);

double deviance_nb (const std::deque<double>& y, const std::deque<double>& mu, const double& phi) {
    double dev=0;
    for (int i=0; i<y.size(); ++i) {
        // We add a small value to protect against zero during division and logging.
        const double& cur_y=y[i]+mildly_low_value;
        const double& cur_mu=mu[i]+mildly_low_value;
        const double product=cur_mu*phi;
        // Calculating the deviance using either the Poisson (small phi*mu), the Gamma (large) or NB (everything else).
        if (product < one_millionth) {
            dev+=cur_y * std::log(cur_y/cur_mu) - (cur_y - cur_mu);
        } else if (product > one_million) {
            dev+=((cur_y - cur_mu)/cur_mu - std::log(cur_y/cur_mu)) * cur_mu/(1+product);
        } else {
            dev+=cur_y * std::log( cur_y/cur_mu ) + (cur_y + 1/phi) * std::log( (cur_mu + 1/phi)/(cur_y + 1/phi) );
        }
    }
    return dev*2;
}

extern "C" {

SEXP mglm_levenberg (SEXP nlib, SEXP ntag, SEXP design, SEXP counts, SEXP disp, SEXP offset, 
        SEXP beta, SEXP fitted, SEXP tol, SEXP maxit){
    const double tolerance=NUMERIC_VALUE(tol);
    const int max_iters=INTEGER_VALUE(maxit);
    PROTECT(design=AS_NUMERIC(design));
    PROTECT(counts=AS_NUMERIC(counts));
    PROTECT(disp=AS_NUMERIC(disp));
    PROTECT(offset=AS_NUMERIC(offset));
    PROTECT(beta=AS_NUMERIC(beta));
    PROTECT(fitted=AS_NUMERIC(fitted));

    // Getting dimensions of the arguments.    
    const long num_tags=(NUMERIC_VALUE(ntag)+0.5);
    const int num_libs=INTEGER_VALUE(nlib);
    const int dlen=LENGTH(design);
    const long clen=LENGTH(counts);
    if (dlen%num_libs!=0) {
        error("Size of design matrix is incompatible with number of libraries.");
    }
    const int num_coefs=dlen/num_libs;
    if (clen!=num_tags*num_libs) { 
        error("Dimensions of the count matrix are not as specified.");
    } else if (LENGTH(beta)!=num_tags*num_coefs) {
        error("Dimensions of the beta matrix do not match to the number of tags and coefficients.");
    } else if (LENGTH(fitted)!=clen) {
        error("Dimensions of the fitted matrix do not match those of the count matrix.");
    } else if (LENGTH(offset)!=clen) {
        error("Dimensions of the offset matrix do not match those of the count matrix.");
    }

    // Initializing pointers.
    std::deque<double*> b_ptrs(num_coefs);
    for (int i=0; i<num_coefs; ++i) { 
        if (i==0) { 
            b_ptrs[i]=NUMERIC_POINTER(beta);
        } else { 
            b_ptrs[i]=b_ptrs[i-1]+num_tags;
        }
    }
    double* design_ptr=NUMERIC_POINTER(design);
    std::deque<double*> c_ptrs(num_libs), f_ptrs(num_libs), o_ptrs(num_libs);
    for (int i=0; i<num_libs; ++i) {
        if (i==0) { 
            c_ptrs[i]=NUMERIC_POINTER(counts); 
            f_ptrs[i]=NUMERIC_POINTER(fitted);
            o_ptrs[i]=NUMERIC_POINTER(offset);
        } else { 
            c_ptrs[i]=c_ptrs[i-1]+num_tags; 
            f_ptrs[i]=f_ptrs[i-1]+num_tags;
            o_ptrs[i]=o_ptrs[i-1]+num_tags;
        }
    }
    double* d_ptr=NUMERIC_POINTER(disp);

    // Initializing output cages. Note that the children of a protected node are also protected.
    SEXP output;
    PROTECT(output=NEW_LIST(5));
    SET_ELEMENT(output, 0, allocMatrix(REALSXP, num_tags, num_coefs));
    SET_ELEMENT(output, 1, allocMatrix(REALSXP, num_tags, num_libs));
    SET_ELEMENT(output, 2, allocVector(REALSXP, num_tags));
    SET_ELEMENT(output, 3, allocVector(INTSXP, num_tags));
    SET_ELEMENT(output, 4, allocVector(LGLSXP, num_tags));

    // Initializing output pointers. Don't have to specify 'AS_NUMERIC' as they are already known.
    std::deque<double*> out_beta_ptrs(num_coefs), out_fitted_ptrs(num_libs);
    for (int i=0; i<num_coefs; ++i) { 
        if (i==0) { out_beta_ptrs[i]=NUMERIC_POINTER(VECTOR_ELT(output, 0)); }
        else { out_beta_ptrs[i]=out_beta_ptrs[i-1]+num_tags; }
    }
    for (int i=0; i<num_libs; ++i) { 
        if (i==0) { out_fitted_ptrs[i]=NUMERIC_POINTER(VECTOR_ELT(output, 1)); }
        else { out_fitted_ptrs[i]=out_fitted_ptrs[i-1]+num_tags; }
    }
    double* dev_ptr=NUMERIC_POINTER(VECTOR_ELT(output, 2));
    int* iter_ptr=INTEGER_POINTER(VECTOR_ELT(output, 3));
    int* fail_ptr=LOGICAL_POINTER(VECTOR_ELT(output, 4));

    // Initializing various arguments.
    std::deque<double> y(num_libs), mu(num_libs), current_beta(num_libs), 
        beta_new(num_libs), mu_new(num_libs);
    double * vx=new double [num_libs*num_coefs];
    double * xvx=new double[num_coefs*num_coefs];
    double * xvx_copy=new double[num_coefs*num_coefs];
    double * dl= new double[num_coefs];
    double * dbeta= new double[num_coefs];
    const char normal='n', transposed='t', uplo='U';
    const double a=1, b=0;
    int info=0, nrhs=1;
    
    // Running through each tag and fitting the NB GLM.
    for (long tag=0; tag<num_tags; ++tag) {
    
        /* Loading the counts and starting means. Note that we do it here
         * to ensure that all the data pointers are incremented; otherwise, 
         * they won't get incremented if they escape into the 'ymax==0'
         * condition, below.
         */
        long ymax=0;
        for (int lib=0; lib<num_libs; ++lib) { 
            const double& current=(y[lib]=*(c_ptrs[lib]++)); 
            mu[lib]=*(f_ptrs[lib]++);
            if (current+0.5>ymax) { ymax=(current+0.5); }
        }  
        for (int coef=0; coef<num_coefs; ++coef) { current_beta[coef]=*(b_ptrs[coef]++); }
        const double& cur_disp=*(d_ptr++);
        double& dev=(*(dev_ptr++)=0);
        int& iter=(*(iter_ptr++)=0);
        int& failed=(*(fail_ptr++)=0);

        // If we start off with all entries at zero, there's really no point continuing. 
        if (ymax==0) {
            for (int coef=0; coef<num_coefs; ++coef) { *(out_beta_ptrs[coef]++)=NA_REAL; }
            for (int lib=0; lib<num_libs; ++lib) { *(out_fitted_ptrs[lib]++)=0; }
            continue;
        }

        // Iterating using reweighted least squares.
		dev=deviance_nb(y, mu, cur_disp);
        double max_info=-1, lambda=0;
        while ((++iter) <= max_iters) {

            // Setting up XVX (or more formally, X*VX).
            for (int row=0; row<num_libs; ++row) {
                const double& cur_mu=mu[row];
                const double cur_prod=cur_mu/(1+cur_mu*cur_disp);
                for (int col=0; col<num_coefs; ++col){ 
                    vx[col*num_libs+row]=design_ptr[col*num_libs+row]*cur_prod;
                }
            }
            F77_NAME(dgemm)(&transposed, &normal, &num_coefs, &num_coefs, &num_libs,
                    &a, design_ptr, &num_libs,
                    vx, &num_libs,
                    &b, xvx, &num_coefs);
            for (int i=0; i<num_coefs; ++i) {
                const double& cur_val=xvx[i*num_coefs+i];
                if (cur_val>max_info) { max_info=cur_val; }
            }
            if (iter==1) {
                lambda=max_info*one_millionth;
                if (lambda < supremely_low_value) { lambda=supremely_low_value; } 
            }

            /* Also setting up 'dl' (derivative of log-likelihoods) as the cross product of the 
             * design matrix with a function of counts/mu/disp. The aim is to repeatedly solve
             * for the step size in a multivariate Newton Raphson using first and second
             * derivatives of the NB likelihood function.
             */
            for (int row=0; row<num_libs; ++row) {
                const double& cur_mu=mu[row];
                double prod=(y[row]-cur_mu)/(1+cur_disp*cur_mu);
                for (int col=0; col<num_coefs; ++col) {
                    if (row==0) { dl[col]=0; }
                    dl[col]+=design_ptr[col*num_libs+row]*prod;
                }
            }

            /* Levenberg/Marquardt damping reduces step size until the deviance increases or no 
             * step can be found that increases the deviance. In short, increases in the deviance
             * are enforced to avoid problems with convergence.
             */ 
            int lev=0;
            bool low_dev=false;
            while (++lev) {
                /* We need to set up copies as the decomposition routine overwrites the originals.
                 * For efficiency, we only refer to the upper triangular for the XVX copy. We also
                 * add 'lambda' to the diagonals (i.e. the bit which stabilises it and helps convergence).
                 * The dl copy is more straightforward. 
                 */
                for (int col=0; col<num_coefs; ++col) {
                    for (int row=0; row<=col; ++row) {
                        const int index=col*num_coefs+row;
                        double& entry=(xvx_copy[index]=xvx[index]);
                        if (row==col) { entry+=lambda; }
                    }
                    dbeta[col]=dl[col];
                }

                // Cholesky decomposition, and then use of the decomposition to solve for Y in (XVX)Y = dl.
                while (1) {                
                    F77_NAME(dpotrf)(&uplo, &num_coefs, xvx_copy, &num_coefs, &info);
                    if (info!=0) {
                        /* If it fails, it MUST mean that the matrix is singular due to numerical imprecision
                         * as all the diagonal entries of the XVX matrix must be positive. This occurs because of 
                         * fitted values being exactly zero; thus, the coefficients attempt to converge to negative 
                         * infinity. This generally forces the step size to be larger (i.e. lambda lower) in order to 
                         * get to infinity faster (which is impossible). Low lambda leads to numerical instability 
                         * and effective singularity. To solve this, we actually increase lambda; this avoids code breakage 
                         * to give the other coefficients a chance to converge. Failure of convergence for the zero-
                         * fitted values isn't a problem as the change in deviance from small --> smaller coefficients isn't 
                         * that great when the true value is negative inifinity.
                         */
                        lambda*=10;
                        for (int col=0; col<num_coefs; ++col) {
                            for (int row=0; row<=col; ++row) {
                                const int index=col*num_coefs+row;
                                double& entry=(xvx_copy[index]=xvx[index]);
                                if (row==col) { entry+=lambda; }
                            }
                        }
                    } else { break; }
                }
                F77_NAME(dpotrs)(&uplo, &num_coefs, &nrhs, xvx_copy, &num_coefs, dbeta, &num_coefs, &info);
                if (info!=0) {
                    std::stringstream err;
                    err << "Solution using Cholesky decomposition failed for tag " << tag+1 << ", check count values.";
                    error(err.str().c_str());
                }

                // Updating beta and the means. 'dbeta' stores 'Y' from the solution of (X*VX)Y=dl, corresponding to a NR step.
                for (int i=0; i<num_coefs; ++i) { beta_new[i]=current_beta[i]+dbeta[i]; }
                for (int lib=0; lib<num_libs; ++lib) {
                    double& cur_mean=(mu_new[lib]=o_ptrs[lib][tag]);
                    for (int coef=0; coef<num_coefs; ++coef) { cur_mean+=design_ptr[coef*num_libs+lib]*beta_new[coef]; }
                    cur_mean=std::exp(cur_mean);
                }

                /* Checking if the deviance has decreased or if it's too small to care about. Either case is good
                 * and means that we'll be using the updated fitted values and coefficients. Otherwise, if we have
                 * to repeat the inner loop, then we want to do so from the original values (as we'll be scaling
                 * lambda up so we want to retake the step from where we were before).
                 */
                const double dev_new=deviance_nb(y, mu_new, cur_disp);
                if (dev_new/ymax < supremely_low_value) { low_dev=true; }
                if (dev_new <= dev || low_dev) {
                    current_beta.swap(beta_new);
                    mu.swap(mu_new);
                    dev=dev_new; 
                    break; 
                }
                
                // Increasing lambda, to increase damping.
                lambda*=2;

                // Excessive damping; steps get so small that it's pointless to continue.
                if (lambda/max_info > 1/supremely_low_value) { 
                    failed=1; 
                    break; 
                }
            } 

            /* Terminating if we failed, if divergence from the exact solution is acceptably low 
             * (cross-product of beta with the log-likelihood derivative) or if the actual deviance 
             * of the fit is acceptably low.
             */
            if (failed) { break; }
            double divergence=0;
            for (int i=0; i<num_coefs; ++i) { divergence+=dl[i]*dbeta[i]; }
            if (low_dev || divergence < tolerance) { break; }

            /* Otherwise, if it is apparent we have not converged i.e. inner loop
             * quits because deviance is decreasing (must be if we survived the breaks),
             * we need bigger steps. To do so, we decrease the damping factor. Note that
             * this only applies if we didn't decrease the damping factor in the inner loop,
             * otherwise we'd just be undoing the work of the inner loop.
             */
            if (lev==1) { lambda/=10; }
        }

        // Storing the data that we've collected, as well as incrementing the output pointers.
        for (int i=0; i<num_coefs; ++i) { *(out_beta_ptrs[i]++)=current_beta[i]; }        
        for (int i=0; i<num_libs; ++i) { *(out_fitted_ptrs[i]++)=mu[i]; }
    }

    delete [] vx;
    delete [] xvx;
    delete [] xvx_copy;
    delete [] dl;
    delete [] dbeta;
    UNPROTECT(7);
    return output;   
}

}

