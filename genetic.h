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

#define MALLOC_FAIL(str) "%s allocation failed!\n", #str

/* e - inverted error condition, f - flag to be returned(not exit) */
#define ASSERT(_e, f, ...)  if (!(_e)) { \
                                fprintf(stderr, __VA_ARGS__); \
                                return f; \
                            }
/* v- variable(name) s-size t-type */
#define ALLOC_ARRAY(v, s, t) v = (t*) malloc(sizeof(t) * s); \
                        ASSERT(v != NULL, FAIL, MALLOC_FAIL(v))
/* v- variable(name) r-rows c-cols t-type */
#define ALLOC_MATRIX(v, r, c, t)  v = (t**) malloc(sizeof(t*) * r); \
                            ASSERT(v != NULL, FAIL, MALLOC_FAIL(v)); \
                            rr=0; \
                            for(;rr < r; ++rr) { \
                                v[rr] = (t*) malloc(sizeof(t) * c); \
                                ASSERT(v[rr] != NULL, FAIL, MALLOC_FAIL(v)); \
                            }

#define _POSIX_C_SOURCE 200809L
#define TRUE 1
#define FALSE 0
#define FAIL -1

typedef enum selection
{
    ROULETTE = 0,
    TOURNAMENT
}selection;

/* TODO: Make an array of ranges */
typedef struct gaconf{
    float *ranges; 
    int dims;       /* dimensions of the target function */
    int size;       /* population size */
    int tour_size;  /* the size of the n/2 tournament */
    float mut_rate; /* probability of mutation from 0 to 1 */
    int elitis;     /* count of units to conserve to elitis */
    int gens;       /* generation count */
    int find_max;   /* 0 if looking for maximum */
    int sel_alg;    /* selection algorithm */
}gaconf;

static int create_matrix(int rows, int cols);
static int create_array(int size);
static void generate_initial_pop(float *ranges, int dims, int size);
static void calculate_costs(float **p, float *c, float(*f)(float*));
static int greater(float a, float b);
static int lesser(float a, float b);
int quicksort_partition(float **p, float *c, int start, int stop);
void quicksort_pop(float **p, float *c, int start, int stop);
static void tournament(int participant_count);
static float prob_max(float cost);
static float prob_min(float cost);
static void roulette(int find_max);
static void crossover_sym();
static void mutation(float mut_rate);
static void elitism(int elitis);
static int ga_init(gaconf *ga);
static float* ga(gaconf* ga, float (*func)(float*));

gaconf def = {NULL, 2, 30, 5, 0.3, 2, 30, FALSE, TOURNAMENT};

int DIMS = 0;
int SIZE = 0;
int MINIMUM = 0;

int (*compare[2]) (float a, float b) = {greater, lesser} ;
void (*sel[2]) (int par) = {roulette, tournament};
float (*prob[2]) (float cost) = {prob_min, prob_max};


static float **pop;
static float *costs;
static float **new_pop;
static float *new_costs;
float *probs;
int *participants;
static int **pairs;

int rr;

static void 
generate_initial_pop(float *ranges, int dims, int size)
{
    int row = 0, col = 0;
    float range, min; 

    for (; row < size; row++){
        col = 0;
        range = ranges[row * 2 + 1] + ranges[row * 2];
        for (; col < dims; col++){
            pop[row][col] = (drand48() * range) + ranges[row*2];
        }
    }
}


static void
calculate_costs(float **p, float *c, float(*f)(float*))
{
    int it = 0;
    
    for(; it < SIZE; it++) {
        c[it] = f(p[it]);
    }
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
quicksort_partition(float **p, float *c, int start, int stop)
{
    int pivot = start;
    int it = start;
    float temp_cost;
    float* temp_chr;
    for (; it < stop; it++) {
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
quicksort_pop(float **p, float *c, int start, int stop)
{
    int pivot; 
    if (start < stop) {
        pivot = quicksort_partition(p, c, start, stop);
        quicksort_pop(p, c, start, pivot-1);
        quicksort_pop(p, c, pivot+1, stop);
    }
}

/* the population is sorted, comparing indices */
static void
tournament(int participant_count)
{
    int pair = 0, it = 0;
    int winners[2] = {INT_MAX, INT_MAX};

    for (; pair < SIZE/2; pair++) {
        participants[0] = drand48() * SIZE;
        it=1;
        for (; it < participant_count; it++) { 
            participants[it] = drand48() * SIZE;
             
            while (participants[it] == participants[it-1]) {
                participants[it] = drand48() * SIZE;
            }
        
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
    
    float prob_sum = 0;

    float range = costs[SIZE-1] - costs[0];

    for (; row < SIZE; row++){
        probs[row] = prob[find_max](costs[row]); 
    }

    for (; pair < SIZE/2; pair++) {
        col = 0;
        for (; col < 2; col++) {
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
}


/* symmetrical crossover */
static void
crossover_sym()
{
    int row = 0, col = 0, pair = 0;
    float r;

    while (row < SIZE) {
        r = (float) drand48();

        for (; col < DIMS; col++) {
            new_pop[row][col] = r * pop[pairs[pair][0]][col];
            new_pop[row][col] += (1-r) * pop[pairs[pair][1]][col];
        }
        
        col = 0;
        row++;

        for (; col < DIMS; col++) {
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
    for (; row < SIZE; row++) {
        col = 0;
        for (; col < DIMS; col++) {  
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
    for (; it < elitis; it++){
        col = 0;
        for (; col < DIMS; col++) {
            new_pop[SIZE-1-it][col] = pop[it][col];
        }
        new_costs[SIZE-1-it] = costs[it];
    }
}


/* dumb function that calls a bunch of mallocs */
static int
ga_init(gaconf *ga)
{
    int pair = 0, it = 0;
 
    SIZE = ga->size;
    DIMS = ga->dims;
    MINIMUM = ga->find_max;

    ALLOC_MATRIX(pop, SIZE, DIMS, float)
    ALLOC_MATRIX(pairs, SIZE/2, 2, int)
    ALLOC_MATRIX(new_pop, SIZE, DIMS, float)

    if(ga->ranges == NULL) {
        ALLOC_ARRAY(ga->ranges, DIMS * 2, float);
        it = 0;
        for (; it < DIMS; it++){
            ga->ranges[it * 2] = -10;
            ga->ranges[it * 2 + 1] = 10;
        }
        printf("Warning: no ranges specified!\n setting as (-10, 10)\n");
    }

    ALLOC_ARRAY(costs, SIZE, float)
    ALLOC_ARRAY(new_costs, SIZE, float)
    
    if (ga->sel_alg == ROULETTE) {
        ALLOC_ARRAY(probs, SIZE, float)
    } else {
        ALLOC_ARRAY(participants, ga->tour_size, int)
    }
    return 0;

}


static float*
ga(gaconf* ga, float (*func)(float*))
{
    float **temp_pop;
    float *temp_costs;
    float prev_best = FLT_MAX;
    ga = ga == NULL ? &def : ga;
    int sel_parameter = ga->sel_alg == ROULETTE ? ga->find_max : ga->tour_size;
    
    ASSERT(ga_init(ga) >= 0, NULL,"ga_init fail!\n") 
    srand48(time(NULL));
    generate_initial_pop(ga->ranges, ga->dims, ga->size);

    calculate_costs(pop, costs, func);
    quicksort_pop(pop, costs, 0, SIZE-1);

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
    }
    
    return pop[0];
}
