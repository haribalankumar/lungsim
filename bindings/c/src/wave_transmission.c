#include "wave_transmission.h"

void evaluate_wave_transmission_c(int *n_time, double *heartrate, double *a0, int *no_freq, double *a, double *b, int *n_bcparams, double *bc_params);

void evaluate_wave_transmission(int n_time, double heartrate, double a0, int no_freq, double *a, double *b, int n_bcparams, double *bc_params)
{
    evaluate_wave_transmission_c(&n_time,&heartrate, &a0, &no_freq, a, b, &n_bcparams,bc_params);
}
