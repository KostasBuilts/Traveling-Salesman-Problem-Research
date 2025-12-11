#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_CITIES 100
#define POPULATION_SIZE 50000
#define GENERATIONS 30
#define MUTATION_RATE 0.01


typedef struct 
{ 
    char name[50]; 
    char comment[50];
    int dimension;
    double x[MAX_CITIES+1];
    double y[MAX_CITIES+1];
} TSPdata;

TSPdata tData;
int num_cities=0;

//dilono tis global metavlites
double cities_distance[MAX_CITIES][MAX_CITIES+1];
int cur_population[POPULATION_SIZE][MAX_CITIES+1];
int new_population[POPULATION_SIZE][MAX_CITIES+1];
int best_route_alltimes[MAX_CITIES+1];
int best_route_population[MAX_CITIES+1];
double fitness[POPULATION_SIZE];
double mean_fitness;
int best_route[21]={0,8,10,6,9,14,12,13,15,11,4,2,1,5,3,20,18,16,19,17,7};





////////////////////////////////////////////////////////////////////////
void strcopy2endline(char *src, char *dst, int start){
    src=src+start;
    while ((*src!=0x0a) && (*src!=0x0d) && (*src!=0x0)){
        *dst++=*src++;
    }
    *dst=0;
}


////////////////////////////////////////////////////////////////////////
void print_city_data()
{
    int i;
    printf("\nNAME: %s", tData.name);
    printf("\nDESC: %s", tData.comment);
    printf("\n");
    printf("+------+----------+----------+\n");
    printf("| city |     X    |     Y    |\n");
    printf("+------+----------+----------+\n");
    for(i=1; i<=tData.dimension; i++){
    
        printf("%6d    %7.2f    %7.2f\n", i, tData.x[i], tData.y[i]);
    }
}




////////////////////////////////////////////////////////////////////////////
int read_filedata(FILE *fptr)
{
    int fHeader=1;
    char line[100], in_buf[50];
    int city=-1;
    double x,y;

    while(fgets(line, 100, fptr) && (city<MAX_CITIES)) {
        if(strstr(line, "NAME:") != NULL) {
            strcopy2endline(line,tData.name,6);
        }

        if(strstr(line, "COMMENT:") != NULL) {
            strcopy2endline(line,tData.comment,9);
        } 

        if(strstr(line, "NODE_COORD_SECTION") != NULL) fHeader=0;

        if (fHeader==0){
            sscanf(line,"%d %lf %lf",  &city, &x, &y);
            tData.x[city]=x;
            tData.y[city]=y;
            tData.dimension=city;
        }
    }
    //printf("read %d of %d cities ", tData.dimension, MAX_CITIES);
    return 1;
}


/////////////////////////////////////////////////////////////////////////
void populate_cites_distance_matrix() 
{
    int i, ii;
    double dx=0, dy=0;

    for (i=1; i<=num_cities; i++) {
        for (ii=1; ii<=num_cities; ii++) {
            dx = tData.x[i] - tData.x[ii];
            dy = tData.y[i] - tData.y[ii];
            cities_distance[i][ii] = sqrt(dx*dx + dy*dy);
        }
    }
}


//////////////////////////////////////////////////////////////////////////
double distance(int a, int b) {
    return cities_distance[a][b];
}


/////////////////////////////////////////////////////////////////////////
double travel_route_length(int *travel_route) {
    double sum = 0;
    int i;
    for (i=1; i<=num_cities-1; i++) {
        sum += distance(travel_route[i], travel_route[i+1]);
    }
    sum += distance(travel_route[num_cities], travel_route[1]);
    return sum;
}


/////////////////////////////////////////////////////////////////////////
void compute_fitness() {
    int i; 
    //calculate the fitness of each individual
    for (i=0; i<POPULATION_SIZE; i++) {
        fitness[i] = 1000.0 / travel_route_length(cur_population[i]);
    }
    //calculate the mean fitness of all population
    for (i=0; i<POPULATION_SIZE; i++) mean_fitness += fitness[i];
    mean_fitness = mean_fitness / POPULATION_SIZE;
}


/////////////////////////////////////////////////////////////////////////
void find_best_fit(int best_route_alltimes[], double *best_fit_alltimes, double *best_fit, int *best_fit_index){
    int i,ii;
    *best_fit=-1;
    *best_fit_index=0;
    for (i=0; i<POPULATION_SIZE; i++) {
        if (fitness[i] > *best_fit) {
            *best_fit = fitness[i];
            *best_fit_index = i;
        }
        if (fitness[i] > *best_fit_alltimes) {
            *best_fit_alltimes = fitness[i];
            for(ii=1; ii<=num_cities; ii++) best_route_alltimes[ii]=cur_population[i][ii];
        } 
    }
}
    

////////////////////////////////////////////////////////////////////////
int select_parent_roulette() {
    double total = 0;
    int i;
    for (i=0; i<POPULATION_SIZE; i++) total += fitness[i];
    double r = ((double) rand() / RAND_MAX) * total;

    double csum = 0;
    for (i=0; i<POPULATION_SIZE; i++) {
        csum += fitness[i];
        if (csum >= r) return i;
    }
    return POPULATION_SIZE - 1;
}


////////////////////////////////////////////////////////////////////////
int select_parent_elitist() {
    static int elites_selected = 0;

    int i;

    // --- 1. Return elites first ---
    if (elites_selected < 2) {
        // Find best individual
        int best = 0;
        for (i=1; i<POPULATION_SIZE; i++) {
            if (fitness[i] > fitness[best]) {
                best = i;
            }
        }
        elites_selected++;
        return best;
    }

    // --- 2. After elites: roulette-wheel selection ---
    double sum = 0.0;
    for (i=0; i<POPULATION_SIZE; i++) {
        sum += fitness[i];
    }

    double r = ((double)rand() / RAND_MAX) * sum;
    double accumulator = 0.0;

    for (i=0; i<POPULATION_SIZE; i++) {
        accumulator += fitness[i];
        if (accumulator >= r) {
            return i;
        }
    }

    return POPULATION_SIZE - 1;
}
////////////////////////////////////////////////////////////////////////////////
void crossover(int *p1, int *p2, int *child) {
    int start, end, i, j, k;
    
    start = rand() % num_cities +1;
    end = rand() % num_cities +1;
    if (start > end) {
        int temp = start;
        start = end;
        end = temp;
    }
    	//init child with all -1
    for (i=1; i<=num_cities; i++) child[i] = -1;
	//copy segment from parent-1
    for (i=start; i<=end; i++) {
        child[i] = p1[i];
    }
	//fill the remain child with parent-2 avoiding dublicates
    k = 1;
    for (i=1; i<=num_cities; i++) {
        int city = p2[i];
        int found = 0;
        for (j=start; j<=end; j++) {
            if (child[j] == city) {
                found = 1;
                break;
            }
        }
        if (!found) {
            while (child[k] != -1) k++;
            child[k] = city;
        }
    }
}


/////////////////////////////////////////////////////////////////////////
void mutate(int *travel_route) {
    double r;
    
    r = (double) rand() / RAND_MAX;
    if (r < MUTATION_RATE) {
        int a = rand() % num_cities +1;
        int b = rand() % num_cities +1;
        int tmp = travel_route[a];
        travel_route[a] = travel_route[b];
        travel_route[b] = tmp;
    }
}


/////////////////////////////////////////////////////////////////////////////
void init_population() {
    int i, j;

    for (i=0; i<POPULATION_SIZE; i++) {
        for (j=1; j<=num_cities; j++) {
            cur_population[i][j] = j;
        }
        for (j=1; j<=num_cities; j++) {
            int a = rand() % num_cities +1;
            int b = rand() % num_cities +1;
            int temp = cur_population[i][a];
            cur_population[i][a] = cur_population[i][b];
            cur_population[i][b] = temp;
        }
    }
}



//////////////////////////////////////////////////////////////////////////
void print_distance_matrix(int start, int end)
{
    int i, ii;

    printf("\n---------\n");
    for (i=start; i<=end; i++) {
        for (ii=start; ii<=end; ii++) {
            printf("%7.2f  ", cities_distance[i][ii]);
        }
        printf("\n");
    }
    printf("-----------\n");
}



//////////////////////////////////////////////////////////////////////////
void print_current_population()
{
    int i, ii;

    printf("\n---------\n");
    for (i=0; i<POPULATION_SIZE; i++) {
        for (ii=1; ii<=num_cities; ii++){
            printf("%d ", cur_population[i][ii]);
        }
        printf("  %f\n", fitness[i]);
    }
    printf("\n---------\n");
}


////////////////////////////////////////////////////////////////////////////
void print_route(int route[])
{
	int i;

    for (i=1; i<=num_cities; i++) {
        printf("%d ", route[i]);
    }

    printf("\nTotal distance: %.3f\n", travel_route_length(route));
}











////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    FILE *fptr;
    int i, ii;
    double tmp;
    double best_fit=-1, best_fit_alltimes=-1;
    int best_index=0;

	//--------------------------------------------------------------------------------------- 
    srand(time(NULL));
    // Open a file in read mode
    fptr = fopen(argv[1], "r");
    if (fptr == NULL) {
            printf ("\n Error: file %s not found\n", argv[1]);
            return(0);
    }
    // read file content and populate tData structure
    if(read_filedata(fptr)==1) printf("\nFile readed succefully\n");
    else printf("\n error during file reading\n");
    fclose(fptr); 
 	//print tData contents
    num_cities=tData.dimension;
    printf("\nNumber of cites = %d\n",num_cities); print_city_data(); 


    //---------------------------------------------------------------------------------------  
    populate_cites_distance_matrix();	// calculate a distance matrix to speed-up executions
    init_population();					// make a random start-up population
    compute_fitness();					// calculate the fitness matrix for each individual
    //print_current_population();
    //print_distance_matrix(0,num_cities);


    //----------------------------------------------------------------------------------------
    printf("\nstarting genetic algorithm...\n");
    printf("-------------------------------------------------------------------------\n");
    printf("  genation    meanfitness   bestfitness-generation   bestfitness-alltimes\n");
    printf("-----------+--------------+------------------------+---------------------\n");
    int gen, pi, idx, k;
    for (gen = 0; gen < GENERATIONS; gen++) {
        for (pi = 0; pi < POPULATION_SIZE; pi++) {
            int p1 = select_parent_elitist();
            int p2 = select_parent_elitist();
            crossover(cur_population[p1], cur_population[p2], new_population[pi]);
            mutate(new_population[pi]);
        }

        for (idx = 0; idx < POPULATION_SIZE; idx++) {
            for (k = 1; k <=num_cities; k++) {
                cur_population[idx][k] = new_population[idx][k];
            }
        }

        compute_fitness();
        find_best_fit(best_route_alltimes, &best_fit_alltimes, &best_fit, &best_index);
        printf(" %5d        %9lf            %9lf              %9lf\n", gen, mean_fitness, best_fit, best_fit_alltimes);     // disply the mean fitness and best fitness
    }

    
    //----------------------------------------------------------------------------------------
    //print_current_population();
    printf("\n-------------------------------------------------");
    printf("\nBest route of last population:\n");
    print_route(cur_population[best_index]);
    printf("\n-------------------------------------------------");
    printf("\nBest route of best generation:\n");
    print_route(best_route_alltimes);  
    printf("\n-------------------------------------------------");
    printf("\nAbsolute Best route:\n");
    print_route(best_route);      

    
    printf("\n\n\n\nobsolute (old stuff)\nBest route found:\n");
    for (i = 0; i < num_cities; i++) {
        printf("%d ", cur_population[best_index][i]);
    }
    printf("\nTotal distance: %.3f\n", travel_route_length(cur_population[best_index]));



    return 0;    

    //1 11 13 14 12 15 9 7 17 19 16 18 10 6 8 5 3 4 2 0  Total distance: 271.947 	
    
    //8 1 5 0 4 2 3 9 7 6 10 18 16 17 19 13 15 12 14 11  Total distance: 269.189   BEST so far !!!!!!!!!!!!!!!!!
    
    
}