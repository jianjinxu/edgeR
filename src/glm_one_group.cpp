#include "glm.h"

std::pair<double,bool> glm_one_group(const int& nlibs, const int& maxit, const double& tolerance, const double* offset, 
#ifdef WEIGHTED		
		const double* weights,	
#endif		
		const double* y, const double* disp, double cur_beta) {
    /* Setting up initial values for beta as the log of the mean of the ratio of counts to offsets.
 	 * This is the exact solution for the gamma distribution (which is the limit of the NB as
 	 * the dispersion goes to infinity. However, if cur_beta is not NA, then we assume it's good. 
 	 */
	bool nonzero=false;
	if (ISNA(cur_beta)) {
		cur_beta=0;
#ifdef WEIGHTED 	   	
 	   	double totweight=0;
#else 
		double totweight=nlibs;
#endif
		for (int j=0; j<nlibs; ++j) {
			const double& cur_val=y[j];
			if (cur_val>low_value) {
#ifdef WEIGHTED			
				cur_beta+=cur_val/std::exp(offset[j]) * weights[j];
#else			
				cur_beta+=cur_val/std::exp(offset[j]);
#endif			
				nonzero=true;
			}
#ifdef WEIGHTED		
			totweight+=weights[j];
#endif			
		}
		cur_beta=std::log(cur_beta/totweight);
	} else {
		for (int j=0; j<nlibs; ++j) { 
			if (y[j] > low_value) { nonzero=true; break; }
		}
	}

	// Skipping to a result for all-zero rows.
	if (!nonzero) { return std::make_pair(R_NegInf, true); }

	// Newton-Raphson iterations to converge to mean.
    bool has_converged=false;
	double dl, info;
	for (int i=0; i<maxit; ++i) {
		dl=0;
 	    info=0;
		for (int j=0; j<nlibs; ++j) {
			const double mu=std::exp(cur_beta+offset[j]), denominator=1+mu*disp[j];
#ifdef WEIGHTED			
			dl+=(y[j]-mu)/denominator * weights[j];
			info+=mu/denominator * weights[j];
#else
			dl+=(y[j]-mu)/denominator;
			info+=mu/denominator;
#endif			
		}
		const double step=dl/info;
		cur_beta+=step;
		if (std::abs(step)<tolerance) {
			has_converged=true;
			break;
		}
	}
	return std::make_pair(cur_beta, has_converged);
}


