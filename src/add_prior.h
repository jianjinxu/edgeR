#ifndef ADD_PRIOR_H
#define ADD_PRIOR_H
#include "matvec_check.h"
#include "utils.h"

class add_prior{
public:
    add_prior(int, int, SEXP, SEXP, bool, bool);
    void fill_and_next();
    const double* const get_priors() const;
    const double* const get_offsets() const;
    const bool same_across_rows() const;
private:
    double* adj_prior, *adj_libs;
    const int num_tags, num_libs;

    matvec_check allp, allo;
    const double* const pptr2, *const optr2;
    const bool logged_in, logged_out;

    int lib, tagdex;
    double ave_lib;
};

#endif
