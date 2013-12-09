#include "utils.h"
extern "C" {
#include "Rmath.h"
}
#include <iostream>

double nbdev (const int& sum, const double& mu, const double& size) {
    const double use_sum=(sum > low_value ? sum : low_value);
   	return 2* ( use_sum*std::log(use_sum/mu) - (use_sum+size)*std::log((use_sum+size)/(mu+size)) );
}

extern "C" {

SEXP R_exact_test_by_deviance(SEXP sums_1, SEXP sums_2, SEXP n_1, SEXP n_2, SEXP disp) try {
	if (!IS_INTEGER(n_1) || LENGTH(n_1)!=1 || !IS_INTEGER(n_2) || LENGTH(n_2)!=1) {
 	   	throw std::runtime_error("number of libraries must be integer scalars"); }
	if (!IS_NUMERIC(disp)) { throw std::runtime_error("dispersion must be a double precision vector"); }
	if (!IS_INTEGER(sums_1) || !IS_INTEGER(sums_2)) { throw std::runtime_error("sums must be integer vectors"); }

    const int n1=INTEGER_VALUE(n_1), n2=INTEGER_VALUE(n_2);
	const int nlibs = n1+n2;
    const int ntags=LENGTH(sums_1);
    if (ntags!=LENGTH(sums_2) || ntags!=LENGTH(disp)) {
        throw std::runtime_error("lengths of input vectors do not match");
    } else if (n1<=0 || n2 <=0) { 
        throw std::runtime_error("number of libraries must be positive for each condition");
    }
    const int* s1_ptr=INTEGER_POINTER(sums_1), *s2_ptr=INTEGER_POINTER(sums_2);
    const double *d_ptr=NUMERIC_POINTER(disp);

    SEXP output=PROTECT(NEW_NUMERIC(ntags));
	try{
		double* p_ptr=NUMERIC_POINTER(output);
    	for (int i=0; i<ntags; ++i) {
        	const int& s1=s1_ptr[i];
 		    const int& s2=s2_ptr[i];
			const int stotal=s1+s2;

			// Computing current means and sizes for each library (probability is the same).
			const double mu = stotal/nlibs;
			const double mu1=mu*n1, mu2=mu*n2, r1=n1/d_ptr[i], r2=n2/d_ptr[i];
			const double p = r1/(r1+mu1);

			/* The aim is to sum conditional probabilities for all partitions of the total sum with deviances 
 			 * greater than that observed for the current partition. We start computing from the extremes
 			 * in both cases.
 			 */
			const double obsdev=nbdev(s1, mu1, r1)+nbdev(s2, mu2, r2);
			double& currentp=(p_ptr[i]=0);
		
			// Going from the left.	
			int j=0;
			while (j <= stotal) {
				if (obsdev <= nbdev(j, mu1, r1)+nbdev(stotal-j, mu2, r2)) { 
					currentp+=dnbinom(j, r1, p, 0) * dnbinom(stotal-j, r2, p, 0);
				} else { break; }
				++j;
			}

			// Going from the right, or what's left of it.
			for (int k=0; k<=stotal-j; ++k) {
				if (obsdev <= nbdev(k, mu2, r2)+nbdev(stotal-k, mu1, r1)) { 
					currentp+=dnbinom(k, r2, p, 0) * dnbinom(stotal-k, r1, p, 0);
				} else { break; }
			}

			const double totalr=r1+r2;
			currentp /= dnbinom(stotal, totalr, totalr/(totalr+mu1+mu2), 0);
		}
	} catch (std::exception& e) { 
		UNPROTECT(1);
		throw;
	}

   	UNPROTECT(1);
    return output;
} catch (std::exception& e) { return mkString(e.what()); }

}
