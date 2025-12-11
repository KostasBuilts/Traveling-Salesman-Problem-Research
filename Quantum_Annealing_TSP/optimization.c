#include "optimization.h"
#include "common.h"
#include <stdlib.h>

// 2-opt local optimization
void two_opt_optimize(TSPInstance *instance, int *tour, double *length) {
    int n = instance->n;
    int improved = 1;
    
    while (improved) {
        improved = 0;
        
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 2; j < n; j++) {
                if (i == 0 && j == n-1) continue;
                
                // Calculate gain
                int a = tour[i];
                int b = tour[(i+1)%n];
                int c = tour[j];
                int d = tour[(j+1)%n];
                
                double gain = instance->dist_matrix[a][b] + instance->dist_matrix[c][d]
                           - instance->dist_matrix[a][c] - instance->dist_matrix[b][d];
                
                if (gain > 1e-9) {
                    // Perform 2-opt swap
                    int start = i + 1;
                    int end = j;
                    while (start < end) {
                        swap_cities(tour, start, end);
                        start++;
                        end--;
                    }
                    *length -= gain;
                    improved = 1;
                }
            }
        }
    }
}

// Validate tour
int validate_tour(int *tour, int n) {
    int *visited = calloc(n, sizeof(int));
    int valid = 1;
    
    for (int i = 0; i < n; i++) {
        if (tour[i] < 0 || tour[i] >= n) {
            valid = 0;
            break;
        }
        visited[tour[i]]++;
    }
    
    if (valid) {
        for (int i = 0; i < n; i++) {
            if (visited[i] != 1) {
                valid = 0;
                break;
            }
        }
    }
    
    free(visited);
    return valid;
}