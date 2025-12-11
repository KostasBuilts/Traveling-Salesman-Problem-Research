#ifndef TSPLIB_H
#define TSPLIB_H

#include "tspinstance.h"

// Parser functions
int parse_tsplib(const char *filename, TSPInstance *instance);
int parse_simple_tsp(const char *filename, TSPInstance *instance);

// Distance calculations
double calculate_distance(TSPInstance *instance, int i, int j);
void build_distance_matrix(TSPInstance *instance);

// Utility functions
double tour_length(TSPInstance *instance, int *tour);
void print_tour(TSPInstance *instance, int *tour);

#endif