#ifndef TSPINSTANCE_H
#define TSPINSTANCE_H

#include "common.h"

// Complete TSPInstance structure
struct TSPInstance 
{
    City *cities;
    int n;
    double **dist_matrix;
    int is_att;
    int is_geo;
    int is_ceiled;
    int is_explicit;
};

// Instance management
void init_tsp_instance(TSPInstance *instance);
void free_tsp_instance(TSPInstance *instance);

#endif