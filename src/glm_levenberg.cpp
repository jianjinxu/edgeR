#include "glm.h"

/* Function to calculate the deviance. Note the protection for very large mu*phi (where we
 * use a gamma instead) or very small mu*phi (where we use the Poisson instead). This 
 * approximation protects against numerical instability introduced by subtrackting
 * a very large log value in (log cur_mu) with another very large logarithm (log cur_mu+1/phi).
 * We need to consider the 'phi' as the approximation is only good when the product is
 * very big or very small.
 */
    
const double one_million=std::pow(10, 6.0), one_millionth=std::pow(10, -6.0);
const double mildly_low_value=std::pow(10, -8.0), supremely_low_value=std::pow(10, -13.0), ridiculously_low_value=std::pow(10, -100.0);

double glm_levenberg::nb_deviance (const double* y, const double* mu, const double& phi) const {
    double dev=0;
    for (int i=0; i<nlibs; ++i) {
        // We add a small value to protect against zero during division and logging.
        const double& cur_y=(y[i] < mildly_low_value ? mildly_low_value : y[i]);
        const double& cur_mu=(mu[i] < mildly_low_value ? mildly_low_value : mu[i]);
        const double product=cur_mu*phi;
        // Calculating the deviance using either the Poisson (small phi*mu), the Gamma (large) or NB (everything else).
        if (product < one_millionth) {
            dev+=cur_y * std::log(cur_y/cur_mu) - (cur_y - cur_mu);
        } else if (product > one_million) {
            dev+=(cur_y - cur_mu)/cur_mu - std::log(cur_y/cur_mu); // * cur_mu/(1+product);
        } else {
            dev+=cur_y * std::log( cur_y/cur_mu ) + (cur_y + 1/phi) * std::log( (cur_mu + 1/phi)/(cur_y + 1/phi) );
        }
    }
    return dev*2;
}

void glm_levenberg::autofill(const double* offset, double* mu, const double* beta) {
	for (int lib=0; lib<nlibs; ++lib) {
		double& cur_mean=(mu[lib]=offset[lib]);
		for (int coef=0; coef<ncoefs; ++coef) { cur_mean+=design[coef*nlibs+lib]*beta[coef]; }
		cur_mean=std::exp(cur_mean);
	}
	return;
}

/* Now, the actual constructors for the GLM object. */

glm_levenberg::glm_levenberg(const int& nl, const int& nc, const double*d, const int& mi, const double& tol) : nlibs(nl), ncoefs(nc),
		maxit(mi), tolerance(tol), info(0) {
	const int len=nlibs*ncoefs;
	design=new double [len];
	for (int i=0; i<len; ++i) { design[i]=d[i]; }

    wx=new double [len];
    xwx=new double [ncoefs*ncoefs];
    xwx_copy=new double [ncoefs*ncoefs];
    dl=new double [ncoefs];
    dbeta=new double [ncoefs];

	mu_new=new double [nlibs];
	beta_new=new double [ncoefs];

	return;
}

glm_levenberg::~glm_levenberg () {
	delete [] design;
	delete [] wx;
	delete [] xwx;
	delete [] xwx_copy;
	delete [] dl;
	delete [] dbeta;
	delete [] mu_new;
	delete [] beta_new;
}

// Now, for the actual fit implementation. 

const char normal='n', transposed='t', uplo='U';
const double a=1, b=0;
const int nrhs=1;

int glm_levenberg::fit(const double* offset, const double* y, const double& disp, double* mu, double* beta) {
	// We expect 'mu' and 'beta' to be supplied. We then check the maximum value of the counts.
    double ymax=0;
    for (int lib=0; lib<nlibs; ++lib) { 
		const double& count=y[lib];
		if (count>ymax) { ymax=count; }
 	}	
    dev=0;
    iter=0;
	failed=false;

    // If we start off with all entries at zero, there's really no point continuing. 
    if (ymax<low_value) {
        for (int coef=0; coef<ncoefs; ++coef) { beta[coef]=NA_REAL; }
        for (int lib=0; lib<nlibs; ++lib) { mu[lib]=0; }
        return 0;
    }
    
	/* Otherwise, we have to make sure 'beta' and 'mu' make sense relative to one another.
 	 * We then proceed to iterating using reweighted least squares.
 	 */
	autofill(offset, mu, beta);
	dev=nb_deviance(y, mu, disp);
    double max_info=-1, lambda=0;

    while ((++iter) <= maxit) {
		for (int i=0; i<ncoefs; ++i) { dl[i]=0; }

		/* Here we set up the matrix XtWX i.e. the Fisher information matrix. X is the design matrix and W is a diagonal matrix
 		 * with the working weights for each observation (i.e. library). The working weights are part of the first derivative of
 		 * the log-likelihood for a given coefficient. When multiplied by two covariates in the design matrix, you get the Fisher 
 		 * information (i.e. variance of the log-likelihood) for that pair. This takes the role of the second derivative of the 
 		 * log-likelihood. The working weights are formed by taking the reciprocal of the product of the variance (in terms of the mean) 
 		 * and the square of the derivative of the link function.
 		 *
 		 * We also set up the actual derivative of the log likelihoods in 'dl'. This is done by multiplying each covariate by the 
 		 * difference between the mu and observation and dividing by the variance and derivative of the link function. This is
 		 * then summed across all observations for each coefficient. The aim is to solve (XtWX)(dbeta)=dl for 'dbeta'. As XtWX
 		 * is the second derivative, and dl is the first, you can see that we are effectively performing a multivariate 
 		 * Newton-Raphson procedure with 'dbeta' as the step.
 		 */
        for (int row=0; row<nlibs; ++row) {
            const double& cur_mu=mu[row];
			const double denom=1+cur_mu*disp;
            const double weight=cur_mu/denom;
			const double deriv=(y[row]-cur_mu)/denom;
            for (int col=0; col<ncoefs; ++col){ 
				const int index=col*nlibs+row;
                wx[index]=design[index]*weight;
				dl[col]+=design[index]*deriv;
            }
        }
        F77_NAME(dgemm)(&transposed, &normal, &ncoefs, &ncoefs, &nlibs,
                &a, design, &nlibs, wx, &nlibs, &b, xwx, &ncoefs);
        for (int i=0; i<ncoefs; ++i) {
            const double& cur_val=xwx[i*ncoefs+i];
            if (cur_val>max_info) { max_info=cur_val; }
        }
        if (iter==1) {
            lambda=max_info*one_millionth;
            if (lambda < supremely_low_value) { lambda=supremely_low_value; } 
        }

        /* Levenberg/Marquardt damping reduces step size until the deviance increases or no 
         * step can be found that increases the deviance. In short, increases in the deviance
         * are enforced to avoid problems with convergence.
         */ 
        int lev=0;
        bool low_dev=false;
        while (++lev) {
			for (int col=0; col<ncoefs; ++col) { dbeta[col]=dl[col]; } // Copying dl to dbeta.
			do {
             	/* We need to set up copies as the decomposition routine overwrites the originals, and 
 				 * we want the originals in case we don't like the latest step. For efficiency, we only 
	 			 * refer to the upper triangular for the XtWX copy (as it should be symmetrical). We also add 
	 			 * 'lambda' to the diagonals. This reduces the step size as the second derivative is increased.
        	     */
         		for (int col=0; col<ncoefs; ++col) {
                	for (int row=0; row<=col; ++row) {
                    	const int index=col*ncoefs+row;
						xwx_copy[index]=xwx[index];
                    	if (row==col) { xwx_copy[index]+=lambda; }
                	}
            	}

            	// Cholesky decomposition, and then use of the decomposition to solve for dbeta in (XtWX)dbeta = dl.
                F77_NAME(dpotrf)(&uplo, &ncoefs, xwx_copy, &ncoefs, &info);
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
                	if (lambda <= 0) { lambda=ridiculously_low_value; } // Just to make sure it actually increases.
                } else { break; }
            } while (1);

            F77_NAME(dpotrs)(&uplo, &ncoefs, &nrhs, xwx_copy, &ncoefs, dbeta, &ncoefs, &info);
            if (info!=0) { return 1; }

            // Updating beta and the means. 'dbeta' stores 'Y' from the solution of (X*VX)Y=dl, corresponding to a NR step.
            for (int i=0; i<ncoefs; ++i) { beta_new[i]=beta[i]+dbeta[i]; }
            autofill(offset, mu_new, beta_new);

            /* Checking if the deviance has decreased or if it's too small to care about. Either case is good
             * and means that we'll be using the updated fitted values and coefficients. Otherwise, if we have
             * to repeat the inner loop, then we want to do so from the original values (as we'll be scaling
             * lambda up so we want to retake the step from where we were before). This is why we don't modify the values
             * in-place until we're sure we want to take the step.
             */
            const double dev_new=nb_deviance(y, mu_new, disp);
            if (dev_new/ymax < supremely_low_value) { low_dev=true; }
            if (dev_new <= dev || low_dev) {
				for (int i=0; i<ncoefs; ++i) { beta[i]=beta_new[i]; }
                for (int i=0; i<nlibs; ++i) { mu[i]=mu_new[i]; }
                dev=dev_new; 
                break; 
            }
            
            // Increasing lambda, to increase damping. Again, we have to make sure it's not zero.
            lambda*=2;
            if (lambda <= 0) { lambda=ridiculously_low_value; }

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
		if (low_dev) { break; }
        double divergence=0;
        for (int i=0; i<ncoefs; ++i) { divergence+=dl[i]*dbeta[i]; }
        if (divergence < tolerance) { break; }

        /* If we quit the inner levenberg loop immediately and survived all the break conditions above, that means that deviance is decreasing
 		 * substantially. Thus, we need larger steps to get there faster. To do so, we decrease the damping factor. Note that this only applies 
 		 * if we didn't decrease the damping factor in the inner levenberg loop, as that would indicate that we need to slow down. 
         */
        if (lev==1) { lambda/=10; }
    }
	return 0;
}

/* Finally, assorted getters. */

const double& glm_levenberg::get_deviance() const {return dev; }

const int& glm_levenberg::get_iterations() const { return iter; }

const bool& glm_levenberg::is_failure() const { return failed; }

