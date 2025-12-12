#include "common.h"
#include "tspinstance.h"
#include "tsplib.h"
#include "quantum.h"
#include "optimization.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int main(int argc, char **argv) {
    // Initialize random seed
    unsigned int seed = time(NULL);
    srand(seed);
    printf("Random seed: %u\n", seed);
    
    // Default OPTIMIZED parameters for quantum annealing
    char *filename = NULL;
    QMCParams params = {
        .classical_temp = 0.1,
        .quantum_temp = 0.5,
        .transverse_field = 5.0,      // Will be auto-adjusted based on problem size
        .num_replicas = 16,           // Will be auto-adjusted
        .num_sweeps = 20,             // Increased for better mixing
        .use_parallel_tempering = 1,  // Use parallel tempering for better exploration
        .num_walkers = 4,             // Will be auto-adjusted
        .num_restarts = 5,            // Multiple restarts for thorough exploration
        .acceptance_target = 0.25     // Target 25% acceptance rate
    };
    
    int use_post_optimization = 1;    // Apply 2-opt after quantum annealing
    int num_threads = omp_get_max_threads();
    int verbose = 1;                  // Print progress information
    double time_limit = 0.0;          // No time limit by default
    int save_tour = 1;                // Save tour to file
    double opt_time = 0.0;            // Initialize opt_time variable
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i+1 < argc) {
            filename = argv[++i];
        }
        else if (strcmp(argv[i], "-temp") == 0 && i+1 < argc) {
            params.classical_temp = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-qtemp") == 0 && i+1 < argc) {
            params.quantum_temp = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-gamma") == 0 && i+1 < argc) {
            params.transverse_field = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-replicas") == 0 && i+1 < argc) {
            params.num_replicas = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-walkers") == 0 && i+1 < argc) {
            params.num_walkers = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-sweeps") == 0 && i+1 < argc) {
            params.num_sweeps = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-restarts") == 0 && i+1 < argc) {
            params.num_restarts = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-threads") == 0 && i+1 < argc) {
            num_threads = atoi(argv[++i]);
            omp_set_num_threads(num_threads);
        }
        else if (strcmp(argv[i], "-time") == 0 && i+1 < argc) {
            time_limit = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-no-pt") == 0) {
            params.use_parallel_tempering = 0;
        }
        else if (strcmp(argv[i], "-no-opt") == 0) {
            use_post_optimization = 0;
        }
        else if (strcmp(argv[i], "-no-save") == 0) {
            save_tour = 0;
        }
        else if (strcmp(argv[i], "-quiet") == 0) {
            verbose = 0;
        }
        else if (strcmp(argv[i], "-seed") == 0 && i+1 < argc) {
            seed = atoi(argv[++i]);
            srand(seed);
            printf("Using specified seed: %u\n", seed);
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Quantum TSP Solver - Enhanced Quantum Annealing\n");
            printf("=================================================\n");
            printf("Usage: %s -f file.tsp [options]\n\n", argv[0]);
            printf("Required:\n");
            printf("  -f file.tsp       TSPLIB instance file\n\n");
            printf("Quantum Annealing Parameters:\n");
            printf("  -temp value       Classical temperature (default: 0.1)\n");
            printf("  -qtemp value      Quantum temperature (default: 0.5)\n");
            printf("  -gamma value      Transverse field strength (default: 5.0)\n");
            printf("  -replicas N       Number of Trotter replicas (default: auto)\n");
            printf("  -walkers N        Number of parallel walkers (default: auto)\n");
            printf("  -sweeps N         Monte Carlo sweeps per step (default: 20)\n");
            printf("  -restarts N       Number of independent restarts (default: 5)\n\n");
            printf("Performance Options:\n");
            printf("  -threads N        Number of CPU threads (default: all=%d)\n", omp_get_max_threads());
            printf("  -time seconds     Maximum runtime in seconds (default: no limit)\n");
            printf("  -seed N           Random seed for reproducibility\n\n");
            printf("Algorithm Options:\n");
            printf("  -no-pt            Disable parallel tempering\n");
            printf("  -no-opt           Disable 2-opt post-optimization\n");
            printf("  -no-save          Don't save tour to file\n");
            printf("  -quiet            Suppress progress output\n\n");
            printf("Information:\n");
            printf("  -h, --help        Show this help message\n");
            printf("\nExamples:\n");
            printf("  %s -f att48.tsp -threads 8 -restarts 10\n", argv[0]);
            printf("  %s -f berlin52.tsp -gamma 10.0 -temp 0.05\n", argv[0]);
            printf("  %s -f test.tsp -no-opt -quiet\n", argv[0]);
            return 0;
        }
        else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            fprintf(stderr, "Use -h for help\n");
            return 1;
        }
    }
    
    if (filename == NULL) {
        fprintf(stderr, "Error: Please specify a TSP file with -f\n");
        fprintf(stderr, "Use -h for help\n");
        return 1;
    }
    
    // Load TSP instance
    TSPInstance instance;
    init_tsp_instance(&instance);
    
    if (verbose) {
        printf("\n=================================================\n");
        printf("Quantum TSP Solver - Enhanced Quantum Annealing\n");
        printf("=================================================\n");
        printf("Loading TSP instance: %s\n", filename);
    }
    
    int success = parse_tsplib(filename, &instance);
    if (!success) {
        if (verbose) printf("Trying simple parser...\n");
        success = parse_simple_tsp(filename, &instance);
    }
    
    if (!success) {
        fprintf(stderr, "Failed to parse TSP file: %s\n", filename);
        fprintf(stderr, "Make sure it's a valid TSPLIB file or a simple coordinate file\n");
        free_tsp_instance(&instance);
        return 1;
    }
    
    if (instance.n <= 1) {
        fprintf(stderr, "Error: Instance must have at least 2 cities\n");
        free_tsp_instance(&instance);
        return 1;
    }
    
    if (verbose) {
        printf("Successfully loaded %d cities\n", instance.n);
        printf("Distance type: ");
        if (instance.is_att) printf("ATT (pseudo-Euclidean)\n");
        else if (instance.is_geo) printf("GEO (geographical)\n");
        else if (instance.is_ceiled) printf("CEIL_2D\n");
        else if (instance.is_explicit) printf("EXPLICIT\n");
        else printf("EUC_2D (standard Euclidean)\n");
    }
    
    // Build distance matrix
    if (verbose) printf("Building distance matrix...\n");
    clock_t build_start = clock();
    build_distance_matrix(&instance);
    clock_t build_end = clock();
    double build_time = (double)(build_end - build_start) / CLOCKS_PER_SEC;
    
    if (verbose) {
        printf("Distance matrix built in %.3f seconds\n", build_time);
        
        // Show some distance samples
        if (instance.n >= 3) {
            printf("Sample distances: ");
            printf("%.2f (0->1), ", instance.dist_matrix[0][1]);
            printf("%.2f (1->2), ", instance.dist_matrix[1][2]);
            if (instance.n > 3) printf("%.2f (2->0)", instance.dist_matrix[2][0]);
            printf("\n");
        }
    }
    
    // Allocate memory for best tour
    int *best_tour = malloc(instance.n * sizeof(int));
    if (!best_tour) {
        fprintf(stderr, "Memory allocation failed for tour\n");
        free_tsp_instance(&instance);
        return 1;
    }
    
    double best_energy;
    
    // Auto-adjust parameters based on problem size if not specified
    if (argc <= 2) { // Only -f was specified, use auto-tuning
        if (verbose) printf("\nAuto-tuning parameters for %d cities...\n", instance.n);
        
        if (instance.n < 10) {
            params.num_replicas = 8;
            params.num_walkers = 2;
            params.num_restarts = 3;
            params.transverse_field = 2.0;
            if (verbose) printf("  Small instance: replicas=%d, walkers=%d, restarts=%d, gamma=%.1f\n",
                              params.num_replicas, params.num_walkers, 
                              params.num_restarts, params.transverse_field);
        }
        else if (instance.n < 30) {
            params.num_replicas = 12;
            params.num_walkers = 4;
            params.num_restarts = 5;
            params.transverse_field = 3.0;
            if (verbose) printf("  Medium instance: replicas=%d, walkers=%d, restarts=%d, gamma=%.1f\n",
                              params.num_replicas, params.num_walkers,
                              params.num_restarts, params.transverse_field);
        }
        else if (instance.n < 100) {
            params.num_replicas = 16;
            params.num_walkers = 6;
            params.num_restarts = 7;
            params.transverse_field = 5.0;
            if (verbose) printf("  Large instance: replicas=%d, walkers=%d, restarts=%d, gamma=%.1f\n",
                              params.num_replicas, params.num_walkers,
                              params.num_restarts, params.transverse_field);
        }
        else {
            params.num_replicas = 20;
            params.num_walkers = 8;
            params.num_restarts = 10;
            params.transverse_field = 8.0;
            if (verbose) printf("  Very large instance: replicas=%d, walkers=%d, restarts=%d, gamma=%.1f\n",
                              params.num_replicas, params.num_walkers,
                              params.num_restarts, params.transverse_field);
        }
    }
    
    if (verbose) {
        printf("\nQuantum Annealing Parameters:\n");
        printf("  Classical temperature: %.3f\n", params.classical_temp);
        printf("  Quantum temperature: %.3f\n", params.quantum_temp);
        printf("  Transverse field (gamma): %.3f\n", params.transverse_field);
        printf("  Replicas (Trotter slices): %d\n", params.num_replicas);
        printf("  Parallel walkers: %d\n", params.num_walkers);
        printf("  Monte Carlo sweeps: %d\n", params.num_sweeps);
        printf("  Independent restarts: %d\n", params.num_restarts);
        printf("  Parallel tempering: %s\n", params.use_parallel_tempering ? "enabled" : "disabled");
        printf("  CPU threads: %d\n", num_threads);
        if (time_limit > 0) printf("  Time limit: %.1f seconds\n", time_limit);
        printf("\nStarting Quantum Annealing...\n");
    }
    
    // Run quantum annealing
    clock_t anneal_start = clock();
    quantum_annealing_optimized(&instance, best_tour, &best_energy, &params);
    clock_t anneal_end = clock();
    
    double anneal_time = (double)(anneal_end - anneal_start) / CLOCKS_PER_SEC;
    
    if (verbose) {
        printf("\nQuantum Annealing completed in %.2f seconds\n", anneal_time);
        printf("Best tour found: %.6f\n", best_energy);
    } else {
        printf("Quantum Annealing: %.2f seconds, length: %.6f\n", anneal_time, best_energy);
    }
    
    // Apply 2-opt local optimization if requested
    if (use_post_optimization) {
        if (verbose) printf("\nApplying 2-opt local optimization...\n");
        
        clock_t opt_start = clock();
        double before_opt = best_energy;
        two_opt_optimize(&instance, best_tour, &best_energy);
        clock_t opt_end = clock();
        
        opt_time = (double)(opt_end - opt_start) / CLOCKS_PER_SEC;
        double improvement = before_opt - best_energy;
        
        if (verbose) {
            printf("2-opt completed in %.3f seconds\n", opt_time);
            printf("Optimized tour length: %.6f\n", best_energy);
            if (improvement > 1e-9) {
                printf("Improvement: %.6f (%.2f%%)\n", improvement, 
                       (improvement / before_opt) * 100.0);
            } else {
                printf("No improvement found (already optimal at this scale)\n");
            }
        }
    }
    
    // Validate the tour
    if (verbose) printf("\nValidating tour...\n");
    
    int valid = validate_tour(best_tour, instance.n);
    
    if (valid) {
        if (verbose) {
            printf("✓ Valid tour found!\n");
            print_tour(&instance, best_tour);
            
            // Calculate exact length to verify
            double exact_length = tour_length(&instance, best_tour);
            printf("Exact tour length: %.6f\n", exact_length);
            
            if (fabs(exact_length - best_energy) > 1e-6) {
                printf("Note: Small rounding difference: %.6f vs %.6f\n", 
                       best_energy, exact_length);
            }
        } else {
            printf("Valid tour: %.6f\n", best_energy);
        }
        
        // Save tour to file if requested
        if (save_tour) {
            char outfile[512];
            snprintf(outfile, sizeof(outfile), "%s.tour", filename);
            FILE *f = fopen(outfile, "w");
            if (f) {
                fprintf(f, "NAME : %s.tour\n", filename);
                fprintf(f, "COMMENT : Found by Enhanced Quantum Annealing\n");
                fprintf(f, "TYPE : TOUR\n");
                fprintf(f, "DIMENSION : %d\n", instance.n);
                fprintf(f, "TOUR_SECTION\n");
                for (int i = 0; i < instance.n; i++) {
                    fprintf(f, "%d\n", best_tour[i] + 1);  // 1-based for TSPLIB
                }
                fprintf(f, "-1\nEOF\n");
                fclose(f);
                
                if (verbose) printf("Tour saved to: %s\n", outfile);
            } else {
                fprintf(stderr, "Warning: Could not save tour to file\n");
            }
        }
        
        // Print summary
        if (verbose) {
            printf("\n=================================================\n");
            printf("SUMMARY\n");
            printf("=================================================\n");
            printf("Instance: %s\n", filename);
            printf("Cities: %d\n", instance.n);
            printf("Total runtime: %.2f seconds\n", build_time + anneal_time + 
                   (use_post_optimization ? opt_time : 0.0));
            printf("Best tour length: %.6f\n", best_energy);
            printf("Random seed: %u\n", seed);
            
            // For known TSPLIB instances, show comparison to known optimal
            if (strstr(filename, "att48") || strstr(filename, "ATT48")) {
                printf("Known optimal for att48: 10628 (difference: %.1f)\n", 
                       best_energy - 10628.0);
            }
            else if (strstr(filename, "berlin52") || strstr(filename, "BERLIN52")) {
                printf("Known optimal for berlin52: 7542 (difference: %.1f)\n", 
                       best_energy - 7542.0);
            }
            else if (strstr(filename, "eil51") || strstr(filename, "EIL51")) {
                printf("Known optimal for eil51: 426 (difference: %.1f)\n", 
                       best_energy - 426.0);
            }
            
            printf("\nTo run again with same random seed:\n");
            printf("  ./quantum_tsp -f %s -seed %u\n", filename, seed);
            printf("=================================================\n");
        }
    } else {
        printf("\n✗ ERROR: Invalid tour generated!\n");
        printf("Tour may have duplicate or out-of-range cities.\n");
        
        // Try to diagnose the problem
        int *visited = calloc(instance.n, sizeof(int));
        if (visited) {
            printf("Tour contents (1-based): ");
            for (int i = 0; i < instance.n && i < 10; i++) {
                printf("%d ", best_tour[i] + 1);
                if (best_tour[i] >= 0 && best_tour[i] < instance.n) {
                    visited[best_tour[i]]++;
                }
            }
            if (instance.n > 10) printf("...");
            printf("\n");
            
            printf("Visit counts (first 10): ");
            for (int i = 0; i < instance.n && i < 10; i++) {
                printf("%d:%d ", i+1, visited[i]);
            }
            printf("\n");
            
            free(visited);
        }
    }
    
    // Cleanup
    free(best_tour);
    free_tsp_instance(&instance);
    
    return valid ? 0 : 1;
}
