#include "utils.h"
#include "glm.h"
extern "C" {
#include "Rmath.h"
}
#ifdef DEBUG
#include <iostream>
#endif

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
			const double phi1=1/r1, phi2=1/r2;
			const double obsdev=compute_unit_nb_deviance(s1, mu1, phi1)+compute_unit_nb_deviance(s2, mu2, phi2);
			double& currentp=(p_ptr[i]=0);
		
			// Going from the left.	
			int j=0;
			while (j <= stotal) {
				if (obsdev <= compute_unit_nb_deviance(j, mu1, phi1)+compute_unit_nb_deviance(stotal-j, mu2, phi2)) { 
					currentp+=dnbinom(j, r1, p, 0) * dnbinom(stotal-j, r2, p, 0);
				} else { break; }
				++j;
			}

			// Going from the right, or what's left of it.
			for (int k=0; k<=stotal-j; ++k) {
				if (obsdev <= compute_unit_nb_deviance(k, mu2, phi2)+compute_unit_nb_deviance(stotal-k, mu1, phi1)) { 
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
