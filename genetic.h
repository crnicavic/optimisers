/*
genetic algorithm library

Copyright (C) 2023 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <limits.h>


#define TRUE 1
#define FALSE 0
#define FAIL -1

typedef enum selection
{
    ROULETTE = 0,
    TOURNAMENT
}selection;

static int greater(float a, float b);
static int lesser(float a, float b);

/* TODO: Make an array of ranges */
typedef struct gaconf{
    float min;      /* the minimum if the range */
    float max;      /* the maximum of the range */
    int dims;       /* dimensions of the target function */
    int size;       /* population size */
    int tour_size;  /* the size of the n/2 tournament */
    float mut_rate; /* probability of mutation from 0 to 1 */
    int elitis;     /* count of units to conserve to elitis */
    float tol;
    int gens;
    int find_max;   /* 0 if looking for maximum */
    int sel_alg;    /* selection algorithm */
}gaconf;


int DIMS = 0;
int SIZE = 0;
int MINIMUM = 0;

int (*compare[2]) (float a, float b);
void (*sel[2]) (int par);
float (*prob[2]) (float cost);

static float **pop;
static float *costs;
static float **new_pop;
static float *new_costs;

int *participants;

static int **pairs;

static void 
generate_initial_pop(float min, float max, int dims, int size)
{
    srand48(time(NULL));
    int row = 0, col = 0;

    pop = (float**) malloc(sizeof(int*) * size);
    if (pop == NULL) {
        printf("Memory allocation failed!\n");
        return;
    }

    for (row = 0; row < size; row++){
        pop[row] = (float*) malloc(sizeof(int) * dims);
        if (pop[row] == NULL) {
            printf("Memory allocation failed!\n");
            return;
        }

        for (col = 0; col < dims; col++){
            pop[row][col] = (drand48() * (max - min)) + min;
        }
    }
}


static void 
print_pop(float** p)
{
    int row = 0, col = 0;
    if (p == NULL || DIMS == 0 || SIZE ==0){
        printf("Can't print empty stuff!\n");
    }
    printf("Population:\n");
    
    for (; row < SIZE; row++) {
        for (col = 0; col < DIMS; col++) {
            printf ("%.2f ", p[row][col]);
        }
        printf("\n");
    }
}


static void
calculate_costs(float** p, float* c, float(*f)(float*))
{
    int it = 0;
    
    for(; it < SIZE; it++) {
        c[it] = f(p[it]);
    }
}


static void
print_costs(float* c)
{
    for (int it = 0; it < SIZE; it++){
        printf("%.2f ", c[it]);
    }
    printf("\n");
}


static int
greater(float a, float b)
{
    return a > b;
}


static int
lesser(float a, float b)
{
    return a < b;
}


/* c - costs, p - it's respective population */
int
quicksort_partition(float** p, float* c, int start, int stop)
{
    /* pivot is practically the place where the pivot element will go */
    int pivot = start;
    int it;
    float temp_cost;
    float* temp_chr;
    for (int it = start; it < stop; it++) {
        if (compare[MINIMUM](c[stop], c[it])) {
            temp_cost = c[it];
            c[it] = c[pivot];
            c[pivot] = temp_cost;
            
            temp_chr = p[it];
            p[it] = p[pivot];
            p[pivot] = temp_chr;

            pivot++;
        }
    }
    temp_cost = c[stop];
    c[stop] = c[pivot];
    c[pivot] = temp_cost;
    
    temp_chr = p[stop];
    p[stop] = p[pivot];
    p[pivot] = temp_chr;

    return pivot;
}


/* it's called quicksort but it sure as shit wont be quick */
void
quicksort_pop(float** p, float* c, int start, int stop)
{
    int pivot; 
    if (start < stop) {
        pivot = quicksort_partition(p, c, start, stop);
        quicksort_pop(p, c, start, pivot-1);
        quicksort_pop(p, c, pivot+1, stop);
    }
}


static void
tournament(int participant_count)
{
    int pair = 0, it = 0;
    int winners[2] = {INT_MAX, INT_MAX};

    if (participants == NULL) {
        printf("Participants memory allocation failed!\n");
        return;
    }

    srand48(time(NULL));

    for (pair = 0; pair < SIZE/2; pair++) {
        participants[0] = drand48() * SIZE;
        for (it = 1; it < participant_count; it++) { 
            participants[it] = drand48() * SIZE;
             
            while (participants[it] == participants[it-1]) {
                participants[it] = drand48() * SIZE;
            }
        
            /* the costs are sorted, so i only need to compare indices */
            if (winners[0] > participants[it]){
                winners[1] = winners[0];
                winners[0] = participants[it];
            } else if (winners[1] > participants[it]){
                winners[1] = participants[it];
            }
        }
        pairs[pair][0] = winners[0]; 
        pairs[pair][1] = winners[1];
    }
}


static float
prob_max(float cost)
{
    return (cost - costs[0]) / (costs[SIZE-1] - costs[0]);
}

static float
prob_min(float cost)
{
    return 1 - (cost - costs[0]) / (costs[SIZE-1] - costs[0]);
}

static void
roulette(int find_max)
{
    float temp;
    int row = 0, col = 0;
    int pair = 0;
    int parent = -1;
    /* TODO: REMOVE THIS MALLOC */
    /* array that holds each probability */
    float* probs = (float*) malloc(sizeof(float) * SIZE);
    
    float prob_sum = 0;

    float range = costs[SIZE-1] - costs[0];

    srand48(time(NULL));

    for (row=0; row < SIZE; row++){
        probs[row] = prob[find_max](costs[row]); 
    }

    for (pair = 0; pair < SIZE/2; pair++) {
        for (col = 0; col < 2; col++) {
            prob_sum = 0;
            temp = drand48();
            parent = 0;
            while(prob_sum < temp){
                prob_sum += probs[parent];
                parent++;
            }
            pairs[pair][col] = parent;
        }
    }
    free(probs);
}


/* symmetrical crossover */
static void
crossover_sym()
{
    int row = 0, col = 0, pair = 0;
    float r;

    while (row < SIZE) {
        r = (float) drand48();

        for ( col = 0; col < DIMS; col++) {
            new_pop[row][col] = r * pop[pairs[pair][0]][col];
            new_pop[row][col] += (1-r) * pop[pairs[pair][1]][col];
        }

        row++;

        for ( col = 0; col < DIMS; col++) {
            new_pop[row][col] = (1-r) * pop[pairs[pair][0]][col]; 
            new_pop[row][col] =+ r * pop[pairs[pair][1]][col];
        }

        row++;
        pair++;
    }
}


static void
mutation(float mut_rate)
{
    float r;
    int row = 0, col = 0;
    for (row = 0; row < SIZE; row++) {
        for (col = 0; col < DIMS; col++) {  
            if (drand48() < mut_rate) {
                r = drand48() * 2;
                new_pop[row][col] += r - 1;
            }
        }
    }
}


static void 
elitism(int elitis)
{
    int it = 0;
    int col = 0;
    for (it = 0; it < elitis; it++){
        for (col = 0; col < DIMS; col++) {
            new_pop[SIZE-1-it][col] = pop[it][col];
        }
        new_costs[SIZE-1-it] = costs[it];
    }
}


/* dumb function that calls a bunch of mallocs */
static int
ga_init(gaconf* ga)
{
    int pair = 0, it = 0;
 
    SIZE = ga->size;
    DIMS = ga->dims;
    MINIMUM = ga->find_max;


    costs = (float*)malloc(sizeof(float) * SIZE);
    if (costs == NULL) {
        printf("Memory allocation failed! [costs]\n");
        return FAIL;
    }
    
    pairs = (int**) malloc(sizeof(int*) * SIZE/2);
    if (pairs == NULL) {
        printf("Memory allocation failed! [pairs]\n");
        return FAIL;
    }

    for ( pair = 0; pair < SIZE/2 ; pair++){
        pairs[pair] = (int*) malloc(sizeof(int) * 2);
        if (pairs[pair] == NULL) {
            printf("Memory allocation failed! [pairs]\n");
            return FAIL;
        }
    }

    new_pop = (float**) malloc(sizeof(float*) * SIZE);
    if (new_pop == NULL) {
        printf("Memory allocation failed! [new_pop]\n");
        return FAIL;
    }

    for (it = 0; it < SIZE; it++) {
        new_pop[it] = (float*) malloc(sizeof(float) * DIMS);
        if (new_pop[it] == NULL) {
            printf("Memory allocation failed! [new_pop]\n");
            return FAIL;
        }
    }

    new_costs = (float*) malloc(sizeof(float) * SIZE);
    if (new_costs == NULL) {
        printf("Memory allocation failed! [new_costs]\n");
        return FAIL;
    }

    participants = (int*) malloc(sizeof(int) * ga->tour_size);
    if (participants == NULL) {
        printf("Memory allocation failed! [participants]\n");
        return FAIL;
    }
    return 0;

}


static float*
ga(gaconf* ga, float (*func)(float*))
{
    float **temp_pop;
    float *temp_costs;
    int sel_parameter = ga->sel_alg == ROULETTE ? ga->find_max : ga->tour_size;
    generate_initial_pop(ga->min, ga->max, ga->dims, ga->size);

    compare[TRUE] = lesser;
    compare[FALSE] = greater;

    sel[ROULETTE] = roulette;
    sel[TOURNAMENT] = tournament;

    prob[TRUE] = prob_max;
    prob[FALSE] = prob_min;

    if (ga_init(ga) < 0) {
        printf("GA INIT FAILED!\n");
        return NULL;
    }

    calculate_costs(pop, costs, func);
    quicksort_pop(pop, costs, 0, SIZE-1);

    printf("%.10f: ", costs[0]); 
    for (int i = 0; i < ga->dims; i++) {
        printf("%.2f ", pop[0][i]);
    }
    puts("");
    for (int it = 0; it < ga->gens; it++) {
        
        sel[ga->sel_alg](sel_parameter);
        crossover_sym();
        mutation(ga->mut_rate);

        calculate_costs(new_pop, new_costs, func);
        quicksort_pop(new_pop, new_costs, 0, SIZE-1); 
        elitism(ga->elitis);
        
        calculate_costs(new_pop, new_costs, func); 
        quicksort_pop(new_pop, new_costs, 0, SIZE-1);

        temp_pop = pop;
        pop = new_pop;
        new_pop = temp_pop;

        temp_costs = costs;
        costs = new_costs;
        new_costs = temp_costs;

        printf("Best: %.10f: ", costs[0]); 
        for (int i = 0; i < ga->dims; i++) {
            printf("%.2f ", pop[0][i]);
        }
        puts("");
    }
    return pop[0];
}

