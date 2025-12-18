#ifndef QUANTUM_H
#define QUANTUM_H

#include "common.h"
#include "tspinstance.h"

typedef struct 
{
    int *tour;
    double energy;
    double *field;
} Replica;

typedef struct 
{
    double classical_temp;
    double quantum_temp;
    double transverse_field;
    int num_replicas;
    int num_sweeps;
    int use_parallel_tempering;
    int num_walkers;
    int num_restarts;
    double acceptance_target;
} QMCParams;

// Quantum annealing with multiple restarts
void quantum_annealing_optimized(TSPInstance *instance, int *best_tour, 
                                double *best_energy, QMCParams *params);

// Wrapper function for backward compatibility
void quantum_annealing(TSPInstance *instance, int *best_tour, double *best_energy,
                      QMCParams *params);

// Adaptive quantum annealing
void adaptive_quantum_annealing(TSPInstance *instance, int *best_tour,
                               double *best_energy, QMCParams *params);

// Replica operations
Replica *create_replica(int n);
void free_replica(Replica *rep);
void reset_replica(Replica *rep, int n);

#endif