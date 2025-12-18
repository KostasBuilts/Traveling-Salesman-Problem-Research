#include "tsplib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <omp.h>

#define PI 3.14159265358979323846

void init_tsp_instance(TSPInstance *instance) 
{
    memset(instance, 0, sizeof(TSPInstance));
}

void free_tsp_instance(TSPInstance *instance) 
{
    if (instance->cities) {
        free(instance->cities);
        instance->cities = NULL;
    }
    if (instance->dist_matrix) {
        for (int i = 0; i < instance->n; i++) {
            free(instance->dist_matrix[i]);
        }
        free(instance->dist_matrix);
        instance->dist_matrix = NULL;
    }
    instance->n = 0;
}

double deg_to_rad(double deg) 
{
    return deg * PI / 180.0;
}

double geo_distance(double lat1, double lon1, double lat2, double lon2) 
{
    double RRR = 6378.388;
    double q1 = cos(lon1 - lon2);
    double q2 = cos(lat1 - lat2);
    double q3 = cos(lat1 + lat2);
    return RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0;
}

double att_distance(double x1, double y1, double x2, double y2) 
{
    double xd = x1 - x2;
    double yd = y1 - y2;
    double rij = sqrt((xd*xd + yd*yd) / 10.0);
    int tij = (int)rij;
    if (tij < rij) return tij + 1.0;
    else return tij;
}

int parse_tsplib(const char *filename, TSPInstance *instance) 
{
    FILE *f = fopen(filename, "r");
    if (!f) 
    {
        fprintf(stderr, "Cannot open file: %s\n", filename);
        return 0;
    }

    char line[2048];
    int dimension = 0;
    int edge_weight_type = 0;
    int found_dimension = 0;
    
    // First pass: read header
    while (fgets(line, sizeof(line), f)) 
    {
        char *ptr = line;
        while (*ptr && isspace(*ptr)) ptr++;
        if (*ptr == '\0' || *ptr == '\n') continue;
        
        char upper_line[2048];
        strncpy(upper_line, line, sizeof(upper_line));
        for (int i = 0; upper_line[i]; i++) 
        {
            upper_line[i] = toupper(upper_line[i]);
        }
        
        if (strstr(upper_line, "DIMENSION")) 
        {
            char *colon = strchr(line, ':');
            if (colon) 
            {
                dimension = atoi(colon + 1);
            } else 
            {
                sscanf(line, "%*s %d", &dimension);
            }
            found_dimension = 1;
        }
        
        if (strstr(upper_line, "EDGE_WEIGHT_TYPE")) 
        {
            if (strstr(upper_line, "ATT")) edge_weight_type = 1;
            else if (strstr(upper_line, "GEO")) edge_weight_type = 2;
            else if (strstr(upper_line, "CEIL_2D")) edge_weight_type = 3;
            else if (strstr(upper_line, "EXPLICIT")) edge_weight_type = 4;
        }
        
        if (strstr(upper_line, "NODE_COORD_SECTION")) 
        {
            break;
        }
    }
    
    if (!found_dimension || dimension <= 0) 
    {
        fclose(f);
        return 0;
    }
    
    instance->n = dimension;
    instance->cities = malloc(dimension * sizeof(City));
    if (!instance->cities) 
    {
        fclose(f);
        return 0;
    }
    
    // Read coordinates
    int cities_read = 0;
    while (fgets(line, sizeof(line), f) && cities_read < dimension) 
    {
        char *ptr = line;
        while (*ptr && isspace(*ptr)) ptr++;
        if (*ptr == '\0' || *ptr == '\n') continue;
        
        if (strstr(line, "EOF")) break;
        
        int idx;
        double x, y;
        if (sscanf(line, "%d %lf %lf", &idx, &x, &y) == 3) 
        {
            if (idx >= 1 && idx <= dimension) 
            {
                instance->cities[idx-1].id = idx;
                instance->cities[idx-1].x = x;
                instance->cities[idx-1].y = y;
                cities_read++;
            }
        }
    }
    
    fclose(f);
    
    if (cities_read != dimension) 
    {
        printf("Warning: Read %d cities, expected %d\n", cities_read, dimension);
        instance->n = cities_read;
    }
    
    instance->is_att = (edge_weight_type == 1);
    instance->is_geo = (edge_weight_type == 2);
    instance->is_ceiled = (edge_weight_type == 3);
    instance->is_explicit = (edge_weight_type == 4);
    
    printf("Loaded %d cities\n", instance->n);
    return 1;
}

int parse_simple_tsp(const char *filename, TSPInstance *instance) 
{
    FILE *f = fopen(filename, "r");
    if (!f) return 0;
    
    // Count lines
    char line[1024];
    int count = 0;
    while (fgets(line, sizeof(line), f)) 
    {
        double x, y;
        if (sscanf(line, "%lf %lf", &x, &y) == 2) 
        {
            count++;
        }
    }
    
    if (count == 0) 
    {
        fclose(f);
        return 0;
    }
    
    rewind(f);
    
    instance->n = count;
    instance->cities = malloc(count * sizeof(City));
    
    count = 0;
    while (fgets(line, sizeof(line), f)) 
    {
        double x, y;
        if (sscanf(line, "%lf %lf", &x, &y) == 2) 
        {
            instance->cities[count].id = count + 1;
            instance->cities[count].x = x;
            instance->cities[count].y = y;
            count++;
        }
    }
    
    fclose(f);
    
    instance->is_att = 0;
    instance->is_geo = 0;
    instance->is_ceiled = 0;
    instance->is_explicit = 0;
    
    printf("Loaded %d cities with simple parser\n", count);
    return 1;
}

double calculate_distance(TSPInstance *instance, int i, int j) 
{
    if (i == j) return 0.0;
    
    City *c1 = &instance->cities[i];
    City *c2 = &instance->cities[j];
    
    if (instance->is_att) 
    {
        return att_distance(c1->x, c1->y, c2->x, c2->y);
    }
    else if (instance->is_geo) 
    {
        double lat1 = deg_to_rad(c1->x);
        double lon1 = deg_to_rad(c1->y);
        double lat2 = deg_to_rad(c2->x);
        double lon2 = deg_to_rad(c2->y);
        return geo_distance(lat1, lon1, lat2, lon2);
    }
    else 
    {
        double dx = c1->x - c2->x;
        double dy = c1->y - c2->y;
        double dist = sqrt(dx*dx + dy*dy);
        if (instance->is_ceiled) 
        {
            return ceil(dist);
        }
        return dist;
    }
}

void build_distance_matrix(TSPInstance *instance) 
{
    int n = instance->n;
    instance->dist_matrix = malloc(n * sizeof(double*));
    
    #pragma omp parallel for
    for (int i = 0; i < n; i++) 
    {
        instance->dist_matrix[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) 
        {
            instance->dist_matrix[i][j] = calculate_distance(instance, i, j);
        }
    }
}

double tour_length(TSPInstance *instance, int *tour) 
{
    double length = 0.0;
    int n = instance->n;
    
    for (int i = 0; i < n; i++) 
    {
        int j = (i + 1) % n;
        length += instance->dist_matrix[tour[i]][tour[j]];
    }
    
    return length;
}

void print_tour(TSPInstance *instance, int *tour) 
{
    printf("Tour (1-based): ");
    for (int i = 0; i < instance->n && i < 20; i++) 
    {
        printf("%d ", tour[i] + 1);
    }
    if (instance->n > 20) printf("...");
    printf("\n");
}