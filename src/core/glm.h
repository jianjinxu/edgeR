#ifndef GLM_H
#define GLM_H

#include "utils.h"

std::pair<double,bool> glm_one_group(const int&, const int&, const double&, const double*, const double*, const double&);

class glm_levenberg {
public:
	glm_levenberg(const int&, const int&, const double*, const int&, const double&);
	~glm_levenberg();
	int fit(const double*, const double*, const double&, double*, double*);

	const bool& is_failure() const;
	const int& get_iterations()  const;
	const double& get_deviance() const;
private:
	const int nlibs;
	const int ncoefs;
	const int maxit;
	const double tolerance;
	
	double* design;
    double * wx;
	double * xwx;
	double * xwx_copy;
	double * dl;
	double * dbeta;
    int info;

	double* mu_new, *beta_new;
	double dev;
	int iter;
	bool failed;

	double nb_deviance(const double*, const double*, const double&) const;
	void autofill(const double*, double*, const double*);
};

class adj_coxreid {
public:
	adj_coxreid(const int&, const int&, const double*);
	~adj_coxreid();
	std::pair<double, bool> compute(const double* wptr);
private:
	const int ncoefs;
	const int nlibs;
	double* design;
	
	double* working_matrix, *work;
	int* pivots;
    int info, lwork;
};

#endif
