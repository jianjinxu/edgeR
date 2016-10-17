#include "utils.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(R_compute_nbdev, 5),
	CALLDEF(R_compute_apl, 6),
	CALLDEF(R_exact_test_by_deviance, 5), 
	CALLDEF(R_loess_by_col, 4),
	CALLDEF(R_maximize_interpolant, 2),

    CALLDEF(R_levenberg, 8),
	CALLDEF(R_get_levenberg_start, 6),
	CALLDEF(R_one_group, 7),
	CALLDEF(R_get_one_way_fitted, 3),
	CALLDEF(R_simple_good_turing, 3),

    CALLDEF(R_add_prior_count, 3),
    CALLDEF(R_calculate_cpm_log, 3),
    CALLDEF(R_calculate_cpm_raw, 2),
    CALLDEF(R_ave_log_cpm, 7),

    CALLDEF(R_check_counts, 1),
    CALLDEF(R_check_finite, 2),
    CALLDEF(R_check_positive, 2),
    CALLDEF(R_check_nonnegative, 2),
    CALLDEF(R_check_zero_fitted, 3),
    CALLDEF(R_check_poisson_bound, 3),
    CALLDEF(R_add_repeat_matrices, 4),
	{NULL, NULL, 0}
};

R_CMethodDef all_c_entries[] = {
    {"processHairpinReads", (DL_FUNC) &processHairpinReads, 23},
    {NULL, NULL, 0}
  };

void attribute_visible R_init_edgeR(DllInfo *dll) {
	R_registerRoutines(dll, all_c_entries, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

}
