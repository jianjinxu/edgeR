#include "add_prior.h"

add_prior::add_prior(int nt, int nl, SEXP priors, SEXP offsets, bool login, bool logout) : num_tags(nt), num_libs(nl),
        allp(priors, num_tags, num_libs), allo(offsets, num_tags, num_libs), pptr2(allp.access()), optr2(allo.access()),
        logged_in(login), logged_out(logout), tagdex(0) {
    adj_prior=(double*)R_alloc(num_libs, sizeof(double));
    adj_libs=(double*)R_alloc(num_libs, sizeof(double));
    return;
}

void add_prior::fill_and_next() {
    if (same_across_rows() && tagdex!=0) {
        // Skipping if all rows are the same, and we've already filled it in once.
        return;
    }

    ave_lib=0;
    for (lib=0; lib<num_libs; ++lib) {
        if (logged_in) { // unlogging to get library sizes, if they were originally logged offsets.
            adj_libs[lib]=std::exp(optr2[lib]);
        } else {
            adj_libs[lib]=optr2[lib];
        }
        ave_lib+=adj_libs[lib];
    }
    ave_lib/=num_libs;

    // Computing the adjusted prior count for each library.
    for (lib=0; lib<num_libs; ++lib) {
        adj_prior[lib]=pptr2[lib]*adj_libs[lib]/ave_lib;
    }

    // Adding it twice back to the library size, and log-transforming.
    for (lib=0; lib<num_libs; ++lib) {
        double& curval=adj_libs[lib];
        curval+=2*adj_prior[lib];
        if (logged_out) {
            curval=std::log(curval);
        }
    }

    ++tagdex;
    allp.advance();
    allo.advance();
    return;
}

const double* const add_prior::get_priors() const { return adj_prior; }

const double* const add_prior::get_offsets()  const { return adj_libs; }

const bool add_prior::same_across_rows() const {
    return (allp.is_row_repeated() && allo.is_row_repeated());
}
