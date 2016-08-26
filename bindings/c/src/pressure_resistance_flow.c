#include "pressure_resistance_flow.h"

void evaluate_prq_c(double *bcinlet, double *bcoutlet, double *targetflow);

void evaluate_prq(double bcinlet, double bcoutlet, double targetflow)
{
    evaluate_prq_c(&bcinlet, &bcoutlet, &targetflow);
}
