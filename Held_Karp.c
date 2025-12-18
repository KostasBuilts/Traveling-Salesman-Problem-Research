#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

/*
 Held–Karp algorithm for the Traveling Salesman Problem (TSP)
 Time complexity: O(n^2 * 2^n)
 Space complexity: O(n * 2^n)
*/

#define MAX_CITIES 20
#define INF 1e18

typedef struct 
{ 
    char name[50]; 
    char comment[50];
    int dimension;
    int cityID[MAX_CITIES];
    double x[MAX_CITIES];
    double y[MAX_CITIES];
} TSPdata;



//dilono tis global metavlites
double dp[1 << MAX_CITIES][MAX_CITIES];
int parent[1 << MAX_CITIES][MAX_CITIES];

TSPdata tData={ 
                "Kiriazis Konstantinos",
                "Internal dataset with 20 cities",
                20, //number of cities
                {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}, // City ID Values
                {10.0, 12.0, 11.0, 13.0,  9.0, 40.0, 41.0, 39.0, 42.0, 38.0, 70.0, 73.0, 75.0, 72.0, 71.0, 20.0, 22.0, 18.0, 21.0, 19.0}, // X Values
                {10.0, 11.0, 14.0,  9.0, 12.0, 40.0, 43.0, 42.0, 39.0, 41.0, 10.0, 11.0,  9.0, 13.0,  8.0, 70.0, 73.0, 69.0, 72.0, 68.0}    // Y Values
            };
int num_cities=0;
double cities_distance[MAX_CITIES][MAX_CITIES];
int tour_path[MAX_CITIES+1];




void strcopy2endline(char *src, char *dst, int start)
{
    src=src+start;
    while ((*src!=0x0a) && (*src!=0x0d) && (*src!=0x0))
    {
        *dst++=*src++;
    }
    *dst=0;
}



int read_filedata(FILE *fptr)
{
    int fHeader=1, eof=0;
    char line[100], in_buf[50];
    int row=0;
    int city=-1;
    double x,y;

    while(fgets(line, 100, fptr) && (row<MAX_CITIES)) 
    {
        if(strstr(line, "NAME:") != NULL) 
        {
            strcopy2endline(line,tData.name,6);
            continue;
        }

        if(strstr(line, "COMMENT:") != NULL) 
        {
            strcopy2endline(line,tData.comment,9);
            continue;
        } 

        if(strstr(line, "NODE_COORD_SECTION") != NULL) 
        {
        	fHeader=0; 
        	continue;
        }
        if (fHeader==0)
        {
            sscanf(line,"%d %lf %lf",  &city, &x, &y);
            tData.x[row]=x;
            tData.y[row]=y;
            tData.cityID[row]=city;
            row++;
            tData.dimension=row;
        }
    }
    return row;
}



void populate_cites_distance_matrix() 
{
    int i, ii;
    double dx=0, dy=0;

    for (i=0; i<num_cities; i++) 
    {
        for (ii=0; ii<num_cities; ii++) 
        {
            dx = tData.x[i] - tData.x[ii];
            dy = tData.y[i] - tData.y[ii];
            cities_distance[i][ii] = sqrt(dx*dx + dy*dy);
        }
    }
}


void print_city_data()
{
    int i;
    printf("\n\n====================\nPrinting Input Daset\n====================");
    printf("\nNAME: %s", tData.name);
    printf("\nDESC: %s", tData.comment);
    printf("\n");
    printf("+------+----------+----------+\n");
    printf("| city |     X    |     Y    |\n");
    printf("+------+----------+----------+\n");
    for(i=0; i<tData.dimension; i++)
    {
    
        printf("%6d    %7.2f    %7.2f\n", tData.cityID[i], tData.x[i], tData.y[i]);
    }
}



void compose_path(int mask, int last) 
{
    static int i;
    if (last == -1) return;
    int prev = parent[mask][last];
    if (prev == -1) 
    {
        return;
    }
    compose_path(mask ^ (1 << last), prev);
    tour_path[i] = last;
    i++;
}





int main(int argc, char **argv) 
{
    int i,j;
    FILE *fptr;


    // Oanoigo to arxeio gia diavasma 
    printf("\n\n\nOpening file ..... ");
    fptr = fopen(argv[1], "r");
    if (fptr == NULL) 
    {
            printf ("\nERROR: file %s not found. Please use a valid file with tsp format.\nSwitching to Internal dataset.\n", argv[1]);
    }
    else 
    {
        // diavazo to arxeio kai genizo tin tData structure
        i=read_filedata(fptr); //pose grammes diavasa
        if(i>=MAX_CITIES) printf("Readed %d cities (the MAXIMUM number) from file %s ",i, argv[1]);
        else printf("Readed %d cities from file %s ",i, argv[1]);
        fclose(fptr); 
    }
    num_cities=tData.dimension;
    print_city_data(); //ektypono ta periexomena tis tData structure


    //---------------------------------------------------------------------------------------  
    populate_cites_distance_matrix();	//proypologizo ton pinaka apostaseon metaxi ton poleon



    // ---------------------------- arxizei o Held-karp -------------------------------
    int u,v,mask;
    int fullMask = 1 << num_cities;

    // arxikopoio ton pinaka katastaseon DP
    for (i = 0; i < fullMask; i++)
    {
        for (j=0; j<num_cities; j++)
        {
            dp[i][j] = INF;
            parent[i][j] = -1;
        }
    }

    // arxiki katastasi afetiria i poli 0
    dp[1][0] = 0;

    // dokimazo ola ta yposynola me tis poleis
    for (mask=1; mask<fullMask; mask++) 
    {
        for (u=0; u<num_cities; u++) 
        {
            if (!(mask & (1 << u))) continue;
            if (dp[mask][u] == INF) continue;

            for (v=0; v<num_cities; v++) 
            {
                if (mask & (1 << v)) continue;
                int nextMask = mask | (1 << v);
                double newCost = dp[mask][u] + cities_distance[u][v];
                if (newCost < dp[nextMask][v]) 
                {
                    dp[nextMask][v] = newCost;
                    parent[nextMask][v] = u;
                }
            }
        }
    }
    // ---------------------------------telos o Held Karp ---------------------------

    // epistrefo stin poli 0
    double tour_cost = INF;
    int lastCity = -1;
    int finalMask = fullMask - 1;

    for (i=1; i<num_cities; i++) 
    {
        double cost = dp[finalMask][i] + cities_distance[i][0];
        if (cost<tour_cost) 
        {
            tour_cost = cost;
            lastCity = i;
        }
    }

    compose_path(finalMask, lastCity);

    printf("\n\n===========\n  RESULTS  \n===========\n");
    printf("Optimal Tour Cost: %.6f\n", tour_cost);
    printf("Optimal Tour Path: ");
    printf("%d ", tData.cityID[0]);
    for(i=0; i<num_cities; i++)
    {
        printf("%d ", tData.cityID[tour_path[i]]);
    }
    printf("\n");
    return 0;
}
