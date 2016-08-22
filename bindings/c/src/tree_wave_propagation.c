#include "tree_wave_propagation.h"

void evaluate_wave_propagation_c(int *n_time, double *a0, int *no_freq, double *a, double *b, int *n_bcparams, double *bc_params);

void evaluate_wave_propagation(int n_time, double a0, int no_freq, double *a, double *b, int n_bcparams, double *bc_params)
{
    evaluate_wave_propagation_c(&n_time, &a0, &no_freq, a, b, &n_bcparams,bc_params);
}

