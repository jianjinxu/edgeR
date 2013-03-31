#include "utils.h"
extern "C" {
#include "Rmath.h"
}

double nbdev (const double& sum, const double& mu, const double& size, const bool& deriv=false) {
    const double& use_sum=(sum > low_value ? sum : low_value);
	if (!deriv) {
    	return use_sum*std::log(use_sum/mu) - (use_sum+size)*std::log((use_sum+size)/(mu+size));
	} else {
		return std::log(use_sum/mu) - std::log((use_sum+size)/(mu+size));
	}
}

extern "C" {

SEXP R_exact_test_by_deviance(SEXP sums_1, SEXP sums_2, SEXP n_1, SEXP n_2, SEXP disp, SEXP big, SEXP tol) try {
    // Setting up the inputs.
    if (!IS_INTEGER(sums_1) || !IS_INTEGER(sums_2)) { throw std::runtime_error("sums must be integer vectors"); }
	if (!IS_NUMERIC(disp)) { throw std::runtime_error("dispersion must be a double precision vector"); }
   
    const int n1=INTEGER_VALUE(n_1), n2=INTEGER_VALUE(n_2);
    const int ntags=LENGTH(sums_1);
    if (ntags!=LENGTH(sums_2) || ntags!=LENGTH(disp)) {
        throw std::runtime_error("lengths of input vectors do not match");
    } else if (n1<=0 || n2 <=0) { 
        throw std::runtime_error("number of libraries must be positive for each condition");
    }
    int* s1_ptr=INTEGER_POINTER(sums_1), *s2_ptr=INTEGER_POINTER(sums_2);
    double *d_ptr=NUMERIC_POINTER(disp);
	const double nr_tolerance=NUMERIC_VALUE(tol);
	const double big_count=NUMERIC_VALUE(big);

	
    // Setting up the outputs.
    SEXP output;
    PROTECT(output=NEW_NUMERIC(ntags));
	double* p_ptr=NUMERIC_POINTER(output);
	try{
    	// Iterating through the tags.
    	const double prop1=n1/double(n1+n2), prop2=n2/double(n1+n2);
    	for (long i=0; i<ntags; ++i) {
        	const double size1=n1/d_ptr[i], size2=n2/d_ptr[i];
        	const int& s1=s1_ptr[i], s2=s2_ptr[i];
        	const int total=s1+s2;
        	const double mu1=total*prop1, mu2=total*prop2;

			if (std::abs(s1-mu1)/s1 < low_value) {
				/* If our count is equal to our mean, then we can just bail and set 
 			 	 * the p-value at 1. It's not going to get any smaller if the deviance is zero. 
 			 	 */
				p_ptr[i]=1;		
		    	continue;	
			} 

			// Sorting out which direction we want to go in.
        	const double threshold_dev=nbdev(s1, mu1, size1) + nbdev(s2, mu2, size2);
			const bool other_is_up=(s1<mu1);
			const double& right_size=(other_is_up ? size2 : size1);
			const double& left_size=(other_is_up ? size1 : size2);
			const double& right_mu=(other_is_up ? mu2 : mu1);
			const double& left_mu=(other_is_up ? mu1 : mu2);

			/* The deviance is a function with one minimum which is monotonic
 		 	 * decreasing upon approach and increasing upon departure. We use 
 		 	 * a Newton-Raphson search to identify the point on the ``other'' 
 		 	 * side which has deviance closest to our observed deviance. 
 		 	 */
			double x=0;
			double step=100;
			
   	 		/* Note the 'minus' in the gradient for the second term; this is because the
 		 	 * differentiation of 'total-x' needs to be reversed to account for the 
 		 	 * negativeness of the 'x'. Also note that the NR search is safe
 		 	 * because the only solution which has a gradient of zero is when the
 		 	 * means are equal to the counts, and that is considered (above).
 		 	 */
			while (std::abs(step) > nr_tolerance) {
				step=(nbdev(x, right_mu, right_size)+nbdev(total-x, left_mu, left_size)-threshold_dev)/
					(nbdev(x, right_mu, right_size, true)-nbdev(total-x, left_mu, left_size, true));
				x-=step;
				if (x > total || x < 0) { throw std::runtime_error("failure during Newton-Raphson procedure"); }
			}
			
        	double& p_out=(p_ptr[i]=0);
			const int& including=(other_is_up ? s1 : s2);

			/* We check if the mu*disp product is large enough for us to use the fact that
 		 	 * the NB distribution is well approximated by the Gamma. This means that the
 		 	 * conditional NB distribution can then be approximated by the Beta distribution.
 		 	 * Note that we only have to check one of them, because mu1*disp2=mu2*disp2=total*disp.
 		 	 */
			if (mu1/size1 > big_count) {
				const double alpha1=mu1/(1+mu1/size1), alpha2=n2/n1*alpha1;
				const double& left_alpha=(other_is_up ? alpha1 : alpha2);
				const double& right_alpha=(other_is_up ? alpha2 : alpha1);
				p_out=pbeta(including/total, left_alpha, right_alpha, 1, 0)
 			   	   +pbeta((x+0.5)/total, right_alpha, left_alpha, 1, 0);
				continue;
			}

        	/* We use lbeta to avoid over/underflow problems resulting from beta.
		 	 * These go away with lbeta because we end up subtracting by the divisor.
		 	 * The price is some loss of precision as the exponent is moved around.
		 	 * However, the number of digits lost is usually small (~3 for
		 	 * the limit of the double datatype). When it gets large, the non-logged
		 	 * version wouldn't be able to handle it anyway, so it's an okay price to pay.
		 	 */
			const double divisor=lbeta(size1, size2);

			/* If the counts are small enough, we iterate. We bascially go through and include
		 	 * our lower partitions on the side of the observed partition and including the observed 
		 	 * partition (hence the +0.5 in the 'including' definition). We do the same for the
		 	 * 'other' side, but we ignore the closest integer for now.
		 	 */
        	double mult=1;
			for (int j=0; j<=including; ++j) {
				p_out+=std::exp(lbeta(j+left_size, total+right_size-j)-divisor)*mult;
            	mult*=(total-j)/(j+1.0);
			}
			mult=1;
			for (int j=0; j<x-0.5; ++j) {
				p_out+=std::exp(lbeta(j+right_size, total+left_size-j)-divisor)*mult;
            	mult*=(total-j)/(j+1.0);
			}
			// We now examine the closest integer. We hold off until this point just to check
			// that it indeed has higher deviance (to protect against NR inaccuracy).
			const double new_x=std::floor(x+0.5);
        	if (nbdev(new_x, right_mu, right_size)+nbdev(total-new_x, left_mu, left_size) > threshold_dev) {
				p_out+=std::exp(lbeta(new_x+right_size, total+left_size-new_x)-divisor)*mult; 		   			
			}
    	}
	} catch (std::exception& e) { 
		UNPROTECT(1);
		throw;
	}

   	UNPROTECT(1);
    return output;
} catch (std::exception& e) { return mkString(e.what()); }

}
