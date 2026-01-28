#include "common.h"
#include <stdlib.h>
#include <time.h>

void swap_cities(int *tour, int i, int j) 
{
    int temp = tour[i];
    tour[i] = tour[j];
    tour[j] = temp;
}

void random_tour(int *tour, int n) 
{
    for (int i = 0; i < n; i++) tour[i] = i;
    
    for (int i = n-1; i > 0; i--) {
        int j = rand() % (i + 1);
        swap_cities(tour, i, j);
    }
}
