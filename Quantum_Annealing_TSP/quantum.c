#include "quantum.h"
#include "tsplib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// Create a new replica
Replica *create_replica(int n) {
    Replica *rep = malloc(sizeof(Replica));
    rep->tour = malloc(n * sizeof(int));
    rep->field = calloc(n, sizeof(double));
    rep->energy = 0.0;
    return rep;
}

void free_replica(Replica *rep) {
    if (rep) {
        free(rep->tour);
        free(rep->field);
        free(rep);
    }
}

void reset_replica(Replica *rep, int n) {
    random_tour(rep->tour, n);
    rep->energy = 0.0;
    for (int i = 0; i < n; i++) {
        rep->field[i] = (rand()/(double)RAND_MAX) > 0.5 ? 1.0 : 0.0;
    }
}

// Different types of moves for better exploration
typedef enum {
    MOVE_SWAP,
    MOVE_2OPT,
    MOVE_3OPT,
    MOVE_REVERSE,
    MOVE_INSERT,
    MOVE_QUANTUM_TUNNEL
} MoveType;

// Try a 3-opt move (more powerful than 2-opt)
int try_3opt_move(TSPInstance *instance, int *tour, int i, int j, int k) {
    int n = instance->n;
    if (i >= j || j >= k) return 0;
    
    // There are 7 possible ways to reconnect 3 broken edges
    // We'll try a few of the most promising ones
    
    int a = tour[i];
    int b = tour[(i+1)%n];
    int c = tour[j];
    int d = tour[(j+1)%n];
    int e = tour[k];
    int f = tour[(k+1)%n];
    
    double current = instance->dist_matrix[a][b] + 
                    instance->dist_matrix[c][d] + 
                    instance->dist_matrix[e][f];
    
    // Try different reconnections
    double best_gain = 0.0;
    int best_type = 0;
    
    // Type 1: Standard 3-opt
    double new1 = instance->dist_matrix[a][c] + 
                 instance->dist_matrix[b][e] + 
                 instance->dist_matrix[d][f];
    double gain1 = current - new1;
    if (gain1 > best_gain) {
        best_gain = gain1;
        best_type = 1;
    }
    
    // Type 2: Another 3-opt variant
    double new2 = instance->dist_matrix[a][d] + 
                 instance->dist_matrix[e][b] + 
                 instance->dist_matrix[c][f];
    double gain2 = current - new2;
    if (gain2 > best_gain) {
        best_gain = gain2;
        best_type = 2;
    }
    
    // Type 3: Yet another variant
    double new3 = instance->dist_matrix[a][e] + 
                 instance->dist_matrix[d][b] + 
                 instance->dist_matrix[c][f];
    double gain3 = current - new3;
    if (gain3 > best_gain) {
        best_gain = gain3;
        best_type = 3;
    }
    
    if (best_gain > 1e-9) {
        // Apply the best move
        if (best_type == 1) {
            // Reverse segment between i+1 and j, and between j+1 and k
            int start1 = i + 1;
            int end1 = j;
            while (start1 < end1) {
                swap_cities(tour, start1, end1);
                start1++;
                end1--;
            }
            
            int start2 = j + 1;
            int end2 = k;
            while (start2 < end2) {
                swap_cities(tour, start2, end2);
                start2++;
                end2--;
            }
        }
        // Add implementations for other types as needed
        return 1;
    }
    
    return 0;
}

// Insert city at different position
int try_insert_move(int *tour, int n) {
    int pos1 = rand() % n;
    int pos2 = rand() % n;
    if (pos1 == pos2) return 0;
    
    int city = tour[pos1];
    
    // Remove from pos1
    if (pos1 < pos2) {
        for (int i = pos1; i < pos2; i++) {
            tour[i] = tour[i+1];
        }
        tour[pos2] = city;
    } else {
        for (int i = pos1; i > pos2; i--) {
            tour[i] = tour[i-1];
        }
        tour[pos2] = city;
    }
    
    return 1;
}

// Enhanced QMC step with multiple move types
double enhanced_qmc_step(TSPInstance *instance, Replica *rep, 
                        double beta, double gamma, int step, int total_steps) {
    int n = instance->n;
    double current_energy = rep->energy;
    double acceptance_rate = 0.0;
    int accepted_moves = 0;
    int attempted_moves = 0;
    
    // Adaptive move probabilities based on annealing progress
    double progress = (double)step / total_steps;
    
    for (int sweep = 0; sweep < 20; sweep++) {
        MoveType move_type;
        double r = (double)rand() / RAND_MAX;
        
        // Adjust move probabilities as annealing progresses
        if (progress < 0.3) {
            // Early: more quantum and large moves
            if (r < 0.3) move_type = MOVE_SWAP;
            else if (r < 0.5) move_type = MOVE_2OPT;
            else if (r < 0.7) move_type = MOVE_3OPT;
            else if (r < 0.85) move_type = MOVE_REVERSE;
            else move_type = MOVE_QUANTUM_TUNNEL;
        } else if (progress < 0.7) {
            // Middle: balanced mix
            if (r < 0.4) move_type = MOVE_SWAP;
            else if (r < 0.6) move_type = MOVE_2OPT;
            else if (r < 0.75) move_type = MOVE_INSERT;
            else if (r < 0.9) move_type = MOVE_REVERSE;
            else move_type = MOVE_QUANTUM_TUNNEL;
        } else {
            // Late: more local refinement
            if (r < 0.6) move_type = MOVE_SWAP;
            else if (r < 0.9) move_type = MOVE_2OPT;
            else move_type = MOVE_INSERT;
        }
        
        attempted_moves++;
        
        // Save current state
        int *old_tour = malloc(n * sizeof(int));
        memcpy(old_tour, rep->tour, n * sizeof(int));
        double old_energy = current_energy;
        double delta_energy = 0.0;
        int move_applied = 0;
        
        switch (move_type) {
            case MOVE_SWAP: {
                int i = rand() % n;
                int j = rand() % n;
                if (i == j) break;
                
                // Calculate energy change
                int prev_i = (i - 1 + n) % n;
                int next_i = (i + 1) % n;
                int prev_j = (j - 1 + n) % n;
                int next_j = (j + 1) % n;
                
                // Remove old edges
                delta_energy -= instance->dist_matrix[rep->tour[prev_i]][rep->tour[i]];
                delta_energy -= instance->dist_matrix[rep->tour[i]][rep->tour[next_i]];
                delta_energy -= instance->dist_matrix[rep->tour[prev_j]][rep->tour[j]];
                delta_energy -= instance->dist_matrix[rep->tour[j]][rep->tour[next_j]];
                
                // Swap
                swap_cities(rep->tour, i, j);
                
                // Add new edges
                delta_energy += instance->dist_matrix[rep->tour[prev_i]][rep->tour[i]];
                delta_energy += instance->dist_matrix[rep->tour[i]][rep->tour[next_i]];
                delta_energy += instance->dist_matrix[rep->tour[prev_j]][rep->tour[j]];
                delta_energy += instance->dist_matrix[rep->tour[j]][rep->tour[next_j]];
                
                move_applied = 1;
                break;
            }
                
            case MOVE_2OPT: {
                int i = rand() % n;
                int j = rand() % n;
                if (abs(i - j) < 2) break;
                if (i > j) { int temp = i; i = j; j = temp; }
                
                int a = rep->tour[i];
                int b = rep->tour[(i+1)%n];
                int c = rep->tour[j];
                int d = rep->tour[(j+1)%n];
                
                double old_dist = instance->dist_matrix[a][b] + instance->dist_matrix[c][d];
                double new_dist = instance->dist_matrix[a][c] + instance->dist_matrix[b][d];
                delta_energy = new_dist - old_dist;
                
                if (delta_energy < 0 || (rand()/(double)RAND_MAX) < exp(-beta * delta_energy)) {
                    int start = i + 1;
                    int end = j;
                    while (start < end) {
                        swap_cities(rep->tour, start, end);
                        start++;
                        end--;
                    }
                    move_applied = 1;
                }
                break;
            }
                
            case MOVE_3OPT: {
                int i = rand() % (n/3);
                int j = rand() % (n/3) + (n/3);
                int k = rand() % (n/3) + (2*n/3);
                
                if (try_3opt_move(instance, rep->tour, i, j, k)) {
                    // Energy will be recalculated
                    move_applied = 1;
                }
                break;
            }
                
            case MOVE_INSERT: {
                if (try_insert_move(rep->tour, n)) {
                    move_applied = 1;
                }
                break;
            }
                
            case MOVE_REVERSE: {
                int start = rand() % n;
                int length = rand() % (n/2) + 1;
                int end = (start + length) % n;
                
                // Reverse the segment
                if (start < end) {
                    while (start < end) {
                        swap_cities(rep->tour, start, end);
                        start++;
                        end--;
                    }
                } else {
                    // Handle wrap-around
                    // Simplified: just do a normal reverse
                    start = rand() % n;
                    end = (start + rand() % (n/2)) % n;
                    if (start > end) { int temp = start; start = end; end = temp; }
                    while (start < end) {
                        swap_cities(rep->tour, start, end);
                        start++;
                        end--;
                    }
                }
                move_applied = 1;
                break;
            }
                
            case MOVE_QUANTUM_TUNNEL: {
                // Quantum-inspired move: large perturbation
                int segment_start = rand() % n;
                int segment_len = rand() % (n/4) + n/8;
                
                // Reverse a random segment
                for (int i = 0; i < segment_len/2; i++) {
                    int idx1 = (segment_start + i) % n;
                    int idx2 = (segment_start + segment_len - i - 1) % n;
                    swap_cities(rep->tour, idx1, idx2);
                }
                move_applied = 1;
                break;
            }
        }
        
        if (move_applied) {
            // Recalculate energy if not already calculated
            if (move_type != MOVE_2OPT && move_type != MOVE_SWAP) {
                current_energy = tour_length(instance, rep->tour);
                delta_energy = current_energy - old_energy;
            } else {
                current_energy += delta_energy;
            }
            
            // Metropolis acceptance with quantum field effect
            double quantum_effect = gamma * (1.0 - progress) * ((double)rand()/RAND_MAX - 0.5);
            double total_delta = delta_energy + quantum_effect;
            
            if (total_delta < 0 || (rand()/(double)RAND_MAX) < exp(-beta * total_delta)) {
                // Accept
                rep->energy = current_energy;
                accepted_moves++;
            } else {
                // Reject: restore old tour
                memcpy(rep->tour, old_tour, n * sizeof(int));
                current_energy = old_energy;
            }
        }
        
        free(old_tour);
    }
    
    acceptance_rate = (double)accepted_moves / attempted_moves;
    return acceptance_rate;
}

// Adaptive parallel tempering with better temperature spacing
void adaptive_parallel_tempering(Replica **replicas, double *temperatures, 
                                int num_replicas, double *acceptance_rates) {
    for (int i = 0; i < num_replicas - 1; i++) {
        double beta1 = 1.0 / temperatures[i];
        double beta2 = 1.0 / temperatures[i+1];
        double delta_beta = beta1 - beta2;
        double delta_energy = replicas[i]->energy - replicas[i+1]->energy;
        
        double swap_prob = exp(delta_beta * delta_energy);
        
        if ((double)rand()/RAND_MAX < swap_prob) {
            // Swap tours
            int *temp_tour = replicas[i]->tour;
            replicas[i]->tour = replicas[i+1]->tour;
            replicas[i+1]->tour = temp_tour;
            
            // Swap energies
            double temp_energy = replicas[i]->energy;
            replicas[i]->energy = replicas[i+1]->energy;
            replicas[i+1]->energy = temp_energy;
            
            // Update acceptance rate
            acceptance_rates[i]++;
        }
    }
}

// Adaptive temperature adjustment - FIXED SYNTAX
void adjust_temperatures(double *temperatures, int num_replicas, 
                        double *acceptance_rates, int steps) {
    double target_rate = 0.25;  // Target 25% acceptance rate
    
    for (int i = 0; i < num_replicas - 1; i++) {
        double current_rate = acceptance_rates[i] / steps;
        
        if (current_rate > target_rate * 1.2) {
            // Too many swaps, increase temperature gap
            temperatures[i] *= 1.05;
            temperatures[i+1] *= 0.95;
        } else if (current_rate < target_rate * 0.8) {
            // Too few swaps, decrease temperature gap
            temperatures[i] *= 0.95;
            temperatures[i+1] *= 1.05;
        }
        
        // Keep temperatures within bounds
        if (temperatures[i] < 0.01) temperatures[i] = 0.01;
        if (temperatures[i] > 100.0) temperatures[i] = 100.0;
        if (temperatures[i+1] < 0.01) temperatures[i+1] = 0.01;
        if (temperatures[i+1] > 100.0) temperatures[i+1] = 100.0;
    }
}

// Main optimized quantum annealing function
void quantum_annealing_optimized(TSPInstance *instance, int *best_tour, 
                                double *best_energy, QMCParams *params) {
    int n = instance->n;
    int P = params->num_replicas;
    int num_walkers = params->num_walkers;
    int num_restarts = params->num_restarts;
    
    printf("Optimized Quantum Annealing\n");
    printf("Cities: %d, Replicas: %d, Walkers: %d, Restarts: %d\n", 
           n, P, num_walkers, num_restarts);
    printf("Threads: %d\n", omp_get_max_threads());
    
    // Track global best across all restarts
    double global_best_energy = 1e100;
    int *global_best_tour = malloc(n * sizeof(int));
    
    // Multiple restarts with different random seeds
    for (int restart = 0; restart < num_restarts; restart++) {
        printf("\n=== Restart %d/%d ===\n", restart + 1, num_restarts);
        
        // Different seed for each restart
        unsigned int seed = time(NULL) + restart * 1000;
        srand(seed);
        
        // Initialize replicas
        Replica ***walkers = malloc(num_walkers * sizeof(Replica**));
        double *walker_best_energies = malloc(num_walkers * sizeof(double));
        int **walker_best_tours = malloc(num_walkers * sizeof(int*));
        
        #pragma omp parallel for
        for (int w = 0; w < num_walkers; w++) {
            walkers[w] = malloc(P * sizeof(Replica*));
            walker_best_tours[w] = malloc(n * sizeof(int));
            walker_best_energies[w] = 1e100;
            
            for (int r = 0; r < P; r++) {
                walkers[w][r] = create_replica(n);
                random_tour(walkers[w][r]->tour, n);
                walkers[w][r]->energy = tour_length(instance, walkers[w][r]->tour);
                
                for (int i = 0; i < n; i++) {
                    walkers[w][r]->field[i] = (rand()/(double)RAND_MAX) > 0.5 ? 1.0 : 0.0;
                }
                
                if (walkers[w][r]->energy < walker_best_energies[w]) {
                    walker_best_energies[w] = walkers[w][r]->energy;
                    memcpy(walker_best_tours[w], walkers[w][r]->tour, n * sizeof(int));
                }
            }
        }
        
        // Adaptive temperature schedule
        double *temperatures = malloc(P * sizeof(double));
        double *acceptance_rates = calloc(P-1, sizeof(double));
        
        // Geometric temperature spacing
        double min_temp = params->quantum_temp;
        double max_temp = min_temp * 10.0;
        for (int r = 0; r < P; r++) {
            temperatures[r] = min_temp * pow(max_temp/min_temp, (double)r/(P-1));
        }
        
        // Adaptive annealing schedule
        int total_steps = 2000;
        double gamma_max = params->transverse_field * (1.0 + restart * 0.2);  // Vary by restart
        double gamma_min = 0.001;
        
        double restart_best_energy = 1e100;
        int *restart_best_tour = malloc(n * sizeof(int));
        
        // Main annealing loop
        for (int step = 0; step < total_steps; step++) {
            double progress = (double)step / total_steps;
            
            // Non-linear annealing schedule
            double gamma = gamma_max * exp(-5.0 * progress) + gamma_min;
            double quantum_temp = params->quantum_temp * (1.0 + 5.0 * (1.0 - progress));
            
            // Process walkers in parallel
            #pragma omp parallel for
            for (int w = 0; w < num_walkers; w++) {
                for (int r = 0; r < P; r++) {
                    double beta = 1.0 / (quantum_temp * sqrt(P));  // Adjusted scaling
                    
                    // Enhanced QMC step (ignore return value to avoid unused variable warning)
                    enhanced_qmc_step(instance, walkers[w][r], 
                                     beta, gamma, step, total_steps);
                    
                    // Update walker best
                    if (walkers[w][r]->energy < walker_best_energies[w]) {
                        walker_best_energies[w] = walkers[w][r]->energy;
                        memcpy(walker_best_tours[w], walkers[w][r]->tour, n * sizeof(int));
                    }
                }
                
                // Adaptive parallel tempering
                if (params->use_parallel_tempering && P > 1 && step % 10 == 0) {
                    adaptive_parallel_tempering(walkers[w], temperatures, P, acceptance_rates);
                }
            }
            
            // Adjust temperatures every 100 steps
            if (step % 100 == 0 && step > 0) {
                adjust_temperatures(temperatures, P, acceptance_rates, 100);
                memset(acceptance_rates, 0, (P-1) * sizeof(double));
            }
            
            // Update restart best
            for (int w = 0; w < num_walkers; w++) {
                if (walker_best_energies[w] < restart_best_energy) {
                    restart_best_energy = walker_best_energies[w];
                    memcpy(restart_best_tour, walker_best_tours[w], n * sizeof(int));
                }
            }
            
            // Progress report
            if (step % 200 == 0) {
                printf("  Step %4d: Gamma=%.4f, Temp=%.4f, Best=%.6f\n",
                       step, gamma, quantum_temp, restart_best_energy);
            }
        }
        
        // Update global best
        if (restart_best_energy < global_best_energy) {
            global_best_energy = restart_best_energy;
            memcpy(global_best_tour, restart_best_tour, n * sizeof(int));
            printf("  New global best: %.6f\n", global_best_energy);
        }
        
        // Cleanup for this restart
        free(restart_best_tour);
        free(temperatures);
        free(acceptance_rates);
        
        for (int w = 0; w < num_walkers; w++) {
            for (int r = 0; r < P; r++) {
                free_replica(walkers[w][r]);
            }
            free(walkers[w]);
            free(walker_best_tours[w]);
        }
        free(walkers);
        free(walker_best_tours);
        free(walker_best_energies);
    }
    
    // Return results
    *best_energy = global_best_energy;
    memcpy(best_tour, global_best_tour, n * sizeof(int));
    
    free(global_best_tour);
}

// Wrapper function with default optimal parameters
void quantum_annealing(TSPInstance *instance, int *best_tour, double *best_energy,
                      QMCParams *params) {
    // Use optimized version with good defaults
    QMCParams optimized_params = *params;
    
    // Auto-adjust parameters based on problem size
    int n = instance->n;
    
    if (n < 20) {
        optimized_params.num_replicas = 8;
        optimized_params.num_walkers = 4;
        optimized_params.num_restarts = 3;
        optimized_params.transverse_field = 3.0;
    } else if (n < 50) {
        optimized_params.num_replicas = 12;
        optimized_params.num_walkers = 6;
        optimized_params.num_restarts = 5;
        optimized_params.transverse_field = 5.0;
    } else if (n < 100) {
        optimized_params.num_replicas = 16;
        optimized_params.num_walkers = 8;
        optimized_params.num_restarts = 7;
        optimized_params.transverse_field = 8.0;
    } else {
        optimized_params.num_replicas = 20;
        optimized_params.num_walkers = 10;
        optimized_params.num_restarts = 10;
        optimized_params.transverse_field = 10.0;
    }
    
    optimized_params.num_sweeps = 20;
    optimized_params.acceptance_target = 0.25;
    
    quantum_annealing_optimized(instance, best_tour, best_energy, &optimized_params);
}

// Adaptive quantum annealing function
void adaptive_quantum_annealing(TSPInstance *instance, int *best_tour,
                               double *best_energy, QMCParams *params) {
    // For now, just call the optimized version
    quantum_annealing_optimized(instance, best_tour, best_energy, params);
}