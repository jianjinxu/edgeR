/* Implements the simple version of the Good-Turing frequency estimator in C++.
 * This is based on the C code written by Geoffrey Sampson from Sussex University.
 * It takes in a vector of observed frequencies and another vector of the same
 * length of frequencies (of observed frequencies). The first vector must be 
 * sorted in ascending order. It also takes a numeric scalar which describes
 * the confidence factor.
 */

#include "edgeR.h"

extern "C" {
	
SEXP simple_good_turing (SEXP obs, SEXP freq, SEXP conf) {
	const double confid_factor=NUMERIC_VALUE(conf);
	PROTECT(obs=AS_NUMERIC(obs));
	PROTECT(freq=AS_NUMERIC(freq));
	const long rows=LENGTH(obs);
	if (rows!=LENGTH(freq)) { error("Length of vectors must match."); }
	double* optr=NUMERIC_POINTER(obs);
    double*	fptr=NUMERIC_POINTER(freq);

	// Prefilling various data structures.
	std::deque<long> obs_int;
	double bigN=0;
	double XYs=0, meanX=0, meanY=0, Xsquares=0;
	std::deque<double> log_obs;
	const long last=rows-1;

	for (long i=0; i<rows; ++i) { 
		const long o=long(optr[i]+0.5);
		obs_int.push_back(o);
		bigN+=o*fptr[i];

		// Computing log data.
		const double& x=(i== 0 ? 0 : optr[i-1]);
		const double logO=std::log(optr[i]),
			  logZ=std::log(2.0*fptr[i]/(i==last ? 2*(optr[i]-x) : optr[i+1]-x));
		log_obs.push_back(logO);
		meanX+=logO;
		meanY+=logZ;
		XYs+=logO*logZ;
		Xsquares+=logO*logO;
	}

	// Finalizing some of the computed variables.
	meanX/=rows;
	meanY/=rows;
	XYs-=meanX*meanY*rows;
	Xsquares-=meanX*meanX*rows;
	const double slope=XYs/Xsquares;
//	std::cout << slope << std::endl;

	// Computing other various bits and pieces.
	const double& PZero = ((rows==0 || obs_int.front()!=1) ? 0 : (fptr[0] / double(bigN)));

	// Setting up the output vector.
	SEXP output=PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(output, 0, ScalarReal(PZero));
	SET_VECTOR_ELT(output, 1, NEW_NUMERIC(rows));
	double* out_ptr=NUMERIC_POINTER(VECTOR_ELT(output, 1));

	// Collecting results.
	double bigNprime=0;
	bool indiffValsSeen=false;
	for (long i=0; i<rows; ++i) {
		const long next_obs=obs_int[i]+1;
		const double y = next_obs*std::exp(slope*(std::log(next_obs)-log_obs[i])); // don't need intercept, cancels out.
//		std::cout << "y for " << i << " is " << y << std::endl;
		if (i==last || obs_int[i+1]!=next_obs) { indiffValsSeen=true; }
		if (!indiffValsSeen) {
			const double& next_n=fptr[i+1];
			const double x = next_obs*next_n/fptr[i];
//			std::cout << "x for " << i << " is " << x << std::endl;
//			std::cout << "test factor is " << confid_factor * x *std::sqrt(1/next_n + 1/fptr[i]) << std::endl;
			if (std::abs(x - y) <= confid_factor * x *std::sqrt(1/next_n + 1/fptr[i])) { // Simplified expression.
				indiffValsSeen=true;
			} else { 
				out_ptr[i]=x; 
//				std::cout << "Using x" << std::endl;
			}
		} 
		if (indiffValsSeen) { 
			out_ptr[i]=y; 
//			std::cout << "Using y" << std::endl;
		}
		bigNprime+=out_ptr[i]*fptr[i];
	}

	// Running through them to compute the remaining bit.
//	std::cout << bigNprime << std::endl;
	const double factor=(1.0-PZero)/bigNprime;
	for (long i=0; i<rows; ++i) { out_ptr[i]*=factor; }

	UNPROTECT(3);
	return output;
}

}
