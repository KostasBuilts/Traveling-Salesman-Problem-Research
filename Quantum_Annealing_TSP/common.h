#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>

#define MAX_CITIES 2000

// Forward declaration for TSPInstance
typedef struct TSPInstance TSPInstance;

// Common structures
typedef struct {
    double x, y;
    int id;
} City;

// Common functions
void swap_cities(int *tour, int i, int j);
void random_tour(int *tour, int n);

#endif