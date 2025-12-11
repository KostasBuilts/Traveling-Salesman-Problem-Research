#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "tspinstance.h"

// Local optimization
void two_opt_optimize(TSPInstance *instance, int *tour, double *length);

// Validation
int validate_tour(int *tour, int n);

#endif