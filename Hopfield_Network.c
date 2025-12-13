#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 20
#define A 1000.0      // Increased: Stronger row constraint
#define B 1000.0      // Increased: Stronger column constraint  
#define C 400.0       // Increased: Stronger distance penalty
#define D 300.0       // Decreased: Weaker self-inhibition
#define DT 0.005      // Smaller step for stability
#define ITER 50000    // More iterations

double x[N][N], u[N][N], dist[N][N], cities[N][2];

void init_cities() {
    for (int i = 0; i < N; i++) {
        cities[i][0] = (double)rand() / RAND_MAX;
        cities[i][1] = (double)rand() / RAND_MAX;
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double dx = cities[i][0] - cities[j][0];
            double dy = cities[i][1] - cities[j][1];
            dist[i][j] = sqrt(dx*dx + dy*dy);
        }
}

int read_tsplib(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) return 0;
    
    char line[256];
    int reading = 0, idx, count = 0;
    double x_coord, y_coord;
    
    while (fgets(line, sizeof(line), f) && count < N) {
        if (strstr(line, "NODE_COORD_SECTION")) reading = 1;
        else if (reading && sscanf(line, "%d %lf %lf", &idx, &x_coord, &y_coord) == 3) {
            if (idx >= 1 && idx <= N) {
                cities[idx-1][0] = x_coord;
                cities[idx-1][1] = y_coord;
                count++;
            }
        }
    }
    fclose(f);
    
    if (count == N) {
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                double dx = cities[i][0] - cities[j][0];
                double dy = cities[i][1] - cities[j][1];
                dist[i][j] = sqrt(dx*dx + dy*dy);
            }
        return 1;
    }
    return 0;
}

void init_network() {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            u[i][j] = -0.2 + ((double)rand()/RAND_MAX)*0.4;  // Wider range
            x[i][j] = 0.5*(1 + tanh(u[i][j]));
        }
}

void update() {
    // Store temporary changes
    double new_u[N][N], new_x[N][N];
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double sum1 = 0, sum2 = 0, sum3 = 0;
            
            // Row constraint: only one city per position
            for (int k = 0; k < N; k++) 
                if (k != j) sum1 += x[i][k];
            
            // Column constraint: only one position per city  
            for (int k = 0; k < N; k++)
                if (k != i) sum2 += x[k][j];
            
            // Distance term: cost of connections
            int jp = (j+1)%N, jm = (j-1+N)%N;
            for (int k = 0; k < N; k++)
                if (k != i)
                    sum3 += dist[i][k] * (x[k][jp] + x[k][jm]);
            
            // Update equation
            double du = -A*sum1 - B*sum2 - C*sum3 - D*(x[i][j] - 0.5);
            new_u[i][j] = u[i][j] + DT*du;
            new_x[i][j] = 0.5*(1 + tanh(new_u[i][j]));
        }
    
    // Apply updates
    memcpy(u, new_u, sizeof(u));
    memcpy(x, new_x, sizeof(x));
}

void decode_solution() {
    int tour[N];
    int used[N] = {0};
    
    // For each position, find best city
    for (int pos = 0; pos < N; pos++) {
        double best = -1;
        int best_city = -1;
        
        for (int city = 0; city < N; city++) {
            if (!used[city] && x[city][pos] > best) {
                best = x[city][pos];
                best_city = city;
            }
        }
        
        // If no unused city found, take the best available
        if (best_city == -1) {
            for (int city = 0; city < N; city++) {
                if (x[city][pos] > best) {
                    best = x[city][pos];
                    best_city = city;
                }
            }
        }
        
        tour[pos] = best_city;
        if (best_city != -1) used[best_city] = 1;
    }
    
    // Fix any duplicates (should be rare with good parameters)
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            if (tour[i] == tour[j]) {
                // Find an unused city
                for (int k = 0; k < N; k++) {
                    if (!used[k]) {
                        tour[j] = k;
                        used[k] = 1;
                        break;
                    }
                }
            }
        }
    }
    
    // Calculate total distance
    double total = 0;
    printf("Tour: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", tour[i] + 1);
        int next = tour[(i+1)%N];
        total += dist[tour[i]][next];
    }
    printf("\nTotal distance: %.4f\n", total);
    
    // Verify it's a valid tour
    int valid = 1;
    int check[N] = {0};
    for (int i = 0; i < N; i++) {
        if (tour[i] < 0 || tour[i] >= N) valid = 0;
        else check[tour[i]]++;
    }
    for (int i = 0; i < N; i++)
        if (check[i] != 1) valid = 0;
    
    if (!valid) printf("Warning: Invalid tour!\n");
}

int main(int argc, char **argv) {
    srand(time(NULL));
    
    if (argc > 1 && read_tsplib(argv[1]))
        printf("Loaded %d cities from %s\n", N, argv[1]);
    else {
        printf("Using random cities\n");
        init_cities();
    }
    
    init_network();
    
    // Run with progress indicator
    printf("Running Hopfield network...\n");
    for (int t = 0; t < ITER; t++) {
        update();
        if (t % 5000 == 0) printf(".");
    }
    printf("\n");
    
    decode_solution();
    
    return 0;
}