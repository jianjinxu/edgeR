#include "glm.h"

std::pair<double,bool> glm_one_group(const int& nlibs, const int& maxit, const double& tolerance, const double* offset, const double* y, const double& disp) {
    /* Setting up initial values for beta as the log of the mean of the ratio of counts to offsets.
 	 * This is the exact solution for the gamma distribution (which is the limit of the NB as
 	 * the dispersion goes to infinity.
 	 */
	bool nonzero=false;
	double cur_beta=0;
	for (int j=0; j<nlibs; ++j) {
		const double& cur_val=y[j];
		if (cur_val>low_value) {
			cur_beta+=cur_val/std::exp(offset[j]);
			nonzero=true;
		}
	}
	if (!nonzero) { return std::make_pair(R_NegInf, true); }

	// If we can't cop out of it, we'll do Newton-Raphson iterations instead.
    bool has_converged=false;
	cur_beta=std::log(cur_beta/nlibs);
	for (int i=0; i<maxit; ++i) {
		double dl=0, info=0;
		for (int j=0; j<nlibs; ++j) {
			const double& cur_val=y[j];
			const double mu=std::exp(cur_beta+offset[j]), denominator=1+mu*disp;
			dl+=(cur_val-mu)/denominator;
			info+=mu/denominator;
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


