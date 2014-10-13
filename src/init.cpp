#include "utils.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(R_compute_nbdev, 3),
	CALLDEF(R_cr_adjust, 3),
	CALLDEF(R_exact_test_by_deviance, 5), 
	CALLDEF(R_levenberg, 11),
	CALLDEF(R_loess_by_col, 4),
	CALLDEF(R_maximize_interpolant, 2),
	CALLDEF(R_one_group, 7),
	CALLDEF(R_simple_good_turing, 3),
	{NULL, NULL, 0}
};

R_CMethodDef all_c_entries[] = {
    {"processHairpinReads", (DL_FUNC) &processHairpinReads, 20},
    {NULL, NULL, 0}
  };

void attribute_visible R_init_edgeR(DllInfo *dll) {
	R_registerRoutines(dll, all_c_entries, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

}
