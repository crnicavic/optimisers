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

/* TODO: I don't need to sort */
/* TODO: Make uneven elitis work, and switch it to percentage */
/* TODO: swap macro */

#include "genetic.h"


static int create_matrix(int rows, int cols);
static int create_array(int size);
static void generate_initial_pop(float *ranges);
static void calculate_costs(float(*f)(float*));
static int greater(float a, float b);
static int lesser(float a, float b);
int partition(int start, int stop);
void quicksort_pop(int start, int stop);
static void tournament(int participant_count);
static float prob_max(float cost);
static float prob_min(float cost);
static void roulette(int find_max);
static void crossover_sym();
static void mutation(float mut_rate);
static void elitism(int elitis);
static int ga_init(gaconf *ga);
float* ga(gaconf* ga, float (*func)(float*));


gaconf def = {NULL, 2, 30, 5, 0.3, 2, 30, FALSE, TOURNAMENT};

int DIMS = 0;
int SIZE = 0;
int MINIMUM = 0;

int (*compare[2]) (float a, float b) = {greater, lesser} ;
void (*sel[2]) (int par) = {roulette, tournament};
float (*prob[2]) (float cost) = {prob_min, prob_max};


float **temp_pop;
float *temp_costs;

float temp_cost;
float* temp_chr;

static float **pop;
static float *costs;
static float **new_pop;
static float *new_costs;
float *probs;
int *participants;
static int **pairs;

int rr;

static void 
generate_initial_pop(float *ranges)
{
    int row = 0, col = 0;
    float range, min; 

    for (; row < SIZE; row++){
        col = 0;
        range = ranges[row * 2 + 1] - ranges[row * 2];
        for (; col < DIMS; col++){
            new_pop[row][col] = (drand48() * range) + ranges[row*2];
        }
    }
}


static void
calculate_costs(float(*f)(float*))
{
    int it = 0;
    
    for(; it < SIZE; it++) {
        new_costs[it] = f(new_pop[it]);
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


int
partition(int start, int stop)
{
    int pivot = start;
    int it = start;
    for (; it < stop; it++) {
        if (compare[MINIMUM](new_costs[stop], new_costs[it])) {
            temp_cost = new_costs[it];
            new_costs[it] = new_costs[pivot];
            new_costs[pivot] = temp_cost;
            
            temp_chr = new_pop[it];
            new_pop[it] = new_pop[pivot];
            new_pop[pivot] = temp_chr;

            pivot++;
        }
    }
    temp_cost = new_costs[stop];
    new_costs[stop] = new_costs[pivot];
    new_costs[pivot] = temp_cost;
    
    temp_chr = new_pop[stop];
    new_pop[stop] = new_pop[pivot];
    new_pop[pivot] = temp_chr;

    return pivot;
}

int
find_kth(int start, int stop, int k)
{
    if (start >= stop) {
        return new_costs[k];
    }

    int pivot = drand48() * (stop-start) + start;
    static int pos;
    static int lsize;

    temp_cost = new_costs[stop];
    new_costs[stop] = new_costs[pivot];
    new_costs[pivot] = temp_cost;
    
    temp_chr = new_pop[stop];
    new_pop[stop] = new_pop[pivot];
    new_pop[pivot] = temp_chr;

    pos = partition(start, stop);
    lsize = pos - start + 1;

    if (lsize == k) {
        return pos;
    } else if (pos > k) {
        return find_kth(start, pos-1, k);
    } else {
        return find_kth(pos+1, stop, k-lsize);
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

static void
tournament(int param) 
{
    return;
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
    /* find elitis best units in pop
     * and set them at the first elitis indices
     */
    temp_pop = pop;
    pop = new_pop;
    new_pop = temp_pop;

    temp_costs = costs;
    costs = new_costs;
    new_costs = temp_costs;
    
    find_kth(0, SIZE-1, elitis); 
    /* find elitis worst units in new_pop
     * and place them at first elitis indices
     * a little hacky but should work
     */

    MINIMUM = FALSE;
    find_kth(0, SIZE-1, elitis);
    
    MINIMUM = TRUE;

    temp_pop = pop;
    pop = new_pop;
    new_pop = temp_pop;

    temp_costs = costs;
    costs = new_costs;
    new_costs = temp_costs;
    
    for(; it < elitis; it++) {
        temp_cost = new_costs[it];
        new_costs[it] = costs[it];
        costs[it] = temp_cost;

        temp_chr = new_pop[it];
        new_pop[it] = pop[it];
        pop[it] = temp_chr;
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


float*
ga(gaconf* ga, float (*func)(float*))
{
    float prev_best = FLT_MAX;
    ga = ga == NULL ? &def : ga;
    int sel_parameter = ga->sel_alg == ROULETTE ? ga->find_max : ga->tour_size;
    
    ASSERT(ga_init(ga) >= 0, NULL,"ga_init fail!\n") 
    srand48(time(NULL));
    
    generate_initial_pop(ga->ranges);
    calculate_costs(func);

    temp_pop = pop;
    pop = new_pop;
    new_pop = temp_pop;

    temp_costs = costs;
    costs = new_costs;
    new_costs = temp_costs;
    
    for (int it = 0; it < ga->gens; it++) {
        
        sel[ga->sel_alg](sel_parameter);
        crossover_sym();
        mutation(ga->mut_rate);

        calculate_costs(func);
        elitism(ga->elitis);
        
        calculate_costs(func); 

        temp_pop = pop;
        pop = new_pop;
        new_pop = temp_pop;

        temp_costs = costs;
        costs = new_costs;
        new_costs = temp_costs;
    }
    
    return pop[0];
}
