#include "utils.h"
#ifndef MATVEC_CHECK_H
#define MATVEC_CHECK_H

class matvec_check {
public:
	matvec_check(SEXP, int, int);
	~matvec_check();
	void advance();
	const double* const access() const;
    const bool is_row_repeated() const;
    const bool is_col_repeated() const;
private:
    int nrow, ncol;
    bool repeat_row, repeat_col;
    
	const double* mycheck;
    int index, libdex, maxdex;
	double* temp;
    const double* toget;
};

class count_holder {
public:
    count_holder(SEXP);
    void fill_and_next(double*);
    int get_ntags() const;
    int get_nlibs() const;

    bool is_data_integer() const;
    const int* get_raw_int() const;
    const double* get_raw_double() const;
private:
    const int* yiptr;
    const double* ydptr;
    int num_tags, num_libs;
    bool is_integer;
    int libdex;
};

#endif
