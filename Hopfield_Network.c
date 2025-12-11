#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 20               // number of cities
#define A 500.0           // penalty constants
#define B 500.0
#define C 200.0
#define D 500.0
#define DT 0.01
#define ITER 30000

double x[N][N];           // neuron outputs
double u[N][N];           // neuron potentials
double dist[N][N];        // distance matrix
double cities[N][2];      // city coordinates

double rand01() { return rand() / (double)RAND_MAX; }

void init_cities() {
    for (int i=0;i<N;i++) {
        cities[i][0] = rand01();
        cities[i][1] = rand01();
    }
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++) {
        double dx = cities[i][0] - cities[j][0];
        double dy = cities[i][1] - cities[j][1];
        dist[i][j] = sqrt(dx*dx + dy*dy);
    }
}

void init_network() {
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++) {
        u[i][j] = rand01()*0.2;
        x[i][j] = 0.5;
    }
}

// Hopfield continuous update
void update() {
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++) {

        double sum1 = 0.0;
        for (int k=0;k<N;k++) if (k!=j) sum1 += x[i][k];

        double sum2 = 0.0;
        for (int k=0;k<N;k++) if (k!=i) sum2 += x[k][j];

        double sum3 = 0.0;
        for (int k=0;k<N;k++) {
            int jp = (j+1) % N;
            int jm = (j-1+N) % N;
            sum3 += dist[i][k] * (x[k][jp] + x[k][jm]);
        }

        double du = -A * sum1 - B * sum2 - C * sum3 - D * (x[i][j] - 0.5);
        u[i][j] += DT * du;
        x[i][j] = 0.5 * (1 + tanh(u[i][j]));
    }
}

void decode_solution() {
    printf("Tour (city indices):\n");
    for (int j=0;j<N;j++) {
        double best = -1.0;
        int best_i = -1;
        for (int i=0;i<N;i++) {
            if (x[i][j] > best) {
                best = x[i][j];
                best_i = i;
            }
        }
        printf("%d ", best_i);
    }
    printf("\n");
}

int main() {
    srand(time(NULL));

    init_cities();
    init_network();

    for (int t=0; t<ITER; t++)
        update();

    decode_solution();
    return 0;
}
