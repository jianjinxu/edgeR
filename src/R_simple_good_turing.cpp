/* Implements the simple version of the Good-Turing frequency estimator in C++.
 * This is based on the C code written by Geoffrey Sampson from Sussex University.
 * It takes in a vector of observed frequencies and another vector of the same
 * length of frequencies (of observed frequencies). The first vector must be 
 * sorted in ascending order. It also takes a numeric scalar which describes
 * the confidence factor.
 */

#include "utils.h"

SEXP R_simple_good_turing (SEXP obs, SEXP freq, SEXP conf) try {
	const double confid_factor=asReal(conf);
	if (!isInteger(obs)) { throw std::runtime_error("observations vector must be integral"); }
	if (!isInteger(freq)) { throw std::runtime_error("frequencies vector must be integral"); }
	const int rows=LENGTH(obs);
	if (rows!=LENGTH(freq)) { throw std::runtime_error("length of vectors must match"); }
	int* optr=INTEGER(obs);
    int* fptr=INTEGER(freq);

	// Prefilling various data structures.
	double bigN=0;
	double XYs=0, meanX=0, meanY=0, Xsquares=0;
	double* log_obs=(double*)R_alloc(rows, sizeof(double));
	const long last=rows-1;

	for (long i=0; i<rows; ++i) { 
		const int& o=optr[i];
		bigN+=o*fptr[i];

		// Computing log data.
		const int& x=(i==0 ? 0 : optr[i-1]);
		const double& logO=(log_obs[i]=std::log(double(optr[i])));
		const double logZ=std::log(2*fptr[i]/double(i==last ? 2*(optr[i]-x) : optr[i+1]-x));
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
	const double& PZero = ((rows==0 || optr[0]!=1) ? 0 : (fptr[0] / double(bigN)));

	// Setting up the output vector.
	SEXP output=PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(output, 0, ScalarReal(PZero));
	SET_VECTOR_ELT(output, 1, allocVector(REALSXP, rows));
	double* out_ptr=REAL(VECTOR_ELT(output, 1));
	try{ 

		// Collecting results.
		double bigNprime=0;
		bool indiffValsSeen=false;
		for (long i=0; i<rows; ++i) {
			const int next_obs=optr[i]+1;
			const double y = next_obs*std::exp(slope*(std::log(double(next_obs))-log_obs[i])); // don't need intercept, cancels out.
	//		std::cout << "y for " << i << " is " << y << std::endl;
			if (i==last || optr[i+1]!=next_obs) { indiffValsSeen=true; }
			if (!indiffValsSeen) {
				const int& next_n=fptr[i+1];
				const double x = next_obs*next_n/double(fptr[i]);
	//			std::cout << "x for " << i << " is " << x << std::endl;
	//			std::cout << "test factor is " << confid_factor * x *std::sqrt(1/next_n + 1/fptr[i]) << std::endl;
				if (std::abs(x - y) <= confid_factor * x *std::sqrt(1.0/next_n + 1.0/fptr[i])) { // Simplified expression.
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
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	} 

	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
