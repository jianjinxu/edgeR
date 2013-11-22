#include "utils.h"
#ifndef MATVEC_CHECK_H
#define MATVEC_CHECK_H

class matvec_check {
public:
	matvec_check(const int, const int, SEXP, const char*);
	void advance();
	const double* const access() const;
private:
	const double* mycheck;
	bool isvec, wasnull;
	int nl;
};

#endif
