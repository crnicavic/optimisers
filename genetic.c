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

#include "genetic.h"

static void generate_initial_pop(float *ranges);
static void calculate_costs(float(*f)(float*));
static int greater(float a, float b);
static int lesser(float a, float b);
static int partition(int start, int stop);
static int find_kth(int start, int stop, int k);
static int* max_min(float *arr);
static void prob_max(int* max_min);
static void prob_min(int* max_min);
static inline int spin();
static void roulette(int find_max);
static int tournament(int participant_count); 
static void brackets(int participant_count);
static void crossover_sym();
static void mutation(float mut_rate);
static void elitism(int elitis);
static int ga_init(gaconf *ga);
static void free_ga(gaconf *ga);
float* ga(gaconf* ga, float (*func)(float*));


gaconf def = {NULL, 2, 30, 5, 0.3, 2, 30, FALSE, ROULETTE};

static int DIMS;
static int SIZE;
static int MAXIMUM;

static int (*compare[2]) (float a, float b) = {lesser, greater} ;
static void (*sel[2]) (int par) = {roulette, brackets};
static void (*prob[2]) (int* max_min) = {prob_min, prob_max};

static float **pop;
static float *costs;
static float **new_pop;
static float *new_costs;
static float *probs;
static int *participants;
static int **pairs;


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


static int
partition(int start, int stop)
{
    int pivot = start;
    int it = start;
    for (; it < stop; it++) {
        if (compare[MAXIMUM](new_costs[stop], new_costs[it])) {
            swap(new_costs[it], new_costs[pivot], float); 
            swap(new_pop[it], new_pop[pivot], float*);
            pivot++;
        }
    }
    
    swap(new_costs[stop], new_costs[pivot], float); 
    swap(new_pop[stop], new_pop[pivot], float*);

    return pivot;
}

static int 
find_kth(int start, int stop, int k)
{
    if (start >= stop) {
        return new_costs[k];
    }

    int pivot = drand48() * (stop-start) + start;
    static int pos;
    static int lsize;

    swap(new_costs[stop], new_costs[pivot], float);
    swap(new_pop[stop], new_pop[pivot], float*);

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


static int* 
max_min(float *arr)
{
    int it = 2;
    static int ret[2] = {0, 0};
    static int min, max;
    if (arr[0] > arr[1]) {
        min = 1;
        max = 0;
    } else {
        min = 0;
        max = 1;
    }

    while (it < SIZE-1) {
        if (arr[it] > arr[it+1]) {
            max = arr[max] < arr[it] ? it : max;
            it++;
            min = arr[min] > arr[it] ? it : min;
        } else {
            min = arr[min] > arr[it] ? it : min;
            it++;
            max = arr[max] < arr[it] ? it : max;
        }
        it++;
    }
    ret[0] = max; 
    ret[1] = min;
    return ret;
}


/* the problem is that the array is no longer sorted */
static void
prob_max(int* max_min)
{
    int row = 0;
    float range = costs[max_min[0]] - costs[max_min[1]];
    float min = costs[max_min[1]];
    for (; row < SIZE; row++) {
        probs[row] = (costs[row] - min) / range;    
    }
}

static void
prob_min(int* max_min)
{
    int row=0;
    float range = costs[max_min[0]] - costs[max_min[1]];
    float min = costs[max_min[1]];
    for (; row < SIZE; row++) {
        probs[row] = 1 - (costs[row] - min) / range;    
    }
}


static inline int
spin()
{
    float temp = drand48();
    float prob_sum = 0;
    int parent = 0;
    while(prob_sum < temp){
        prob_sum += probs[parent]; 
        parent++;
    }

    return parent;
}


static void
roulette(int find_max)
{
    int row = 0, col = 0;
    int pair = 0;
    int *r = max_min(costs); 

    prob[find_max](r);

    for (; pair < SIZE/2; pair++) {
        col = 0;
        for (; col < 2; col++) {
            pairs[pair][col] = spin();
        }
    }
}


static int
tournament(int participant_count) 
{
    int it = 1, win = -1;

    participants[0] = drand48() * SIZE;
    win = participants[0];
    for (;it < participant_count; it++) {
        participants[it] = drand48() * SIZE;
        
        while (participants[it] == participants[it-1]) {
            participants[it] = drand48() * SIZE;
        }

        if(compare[MAXIMUM](costs[participants[it]],costs[win])) {
            win = participants[it];
        }
    }
    return win;
}


static void
brackets(int participant_count)
{
    int row = 0, col = 0;
    int pair = 0;
    for (; pair < SIZE/2; pair++) {
        col = 0;
        for (; col < 2; col++) {
            pairs[pair][col] = tournament(participant_count);
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
            r = drand48() < mut_rate ? (drand48() * 2) - 1 : 0;
            new_pop[row][col] += r;
        }
    }
}


static void 
elitism(int elitis)
{
    int it = 0;
   
    swap(pop, new_pop, float**);
    swap(costs, new_costs, float*);

    /* looks for elitis best */
    find_kth(0, SIZE-1, elitis); 

    MAXIMUM = MAXIMUM ^ 1;
    
    /*looks for elitis worst */
    find_kth(0, SIZE-1, elitis);
    
    MAXIMUM = MAXIMUM ^ 1;

    swap(pop, new_pop, float**);
    swap(costs, new_costs, float*);
    
    for(; it < elitis; it++) {
        swap(costs[it], new_costs[it], float);
        swap(pop[it], new_pop[it], float*);
    }
}

/* dumb function that calls a bunch of mallocs */
static int
ga_init(gaconf *ga)
{
    int pair = 0, it = 0;
 
    SIZE = ga->size;
    DIMS = ga->dims;
    MAXIMUM = ga->find_max;

    ALLOC_MATRIX(pop, SIZE, DIMS, float);
    ALLOC_MATRIX(pairs, SIZE/2, 2, int);
    ALLOC_MATRIX(new_pop, SIZE, DIMS, float);

    if(ga->ranges == NULL) {
        ALLOC_ARRAY(ga->ranges, DIMS * 2, float);
        it = 0;
        for (; it < DIMS; it++){
            ga->ranges[it * 2] = -10;
            ga->ranges[it * 2 + 1] = 10;
        }
        printf("Warning: no ranges specified!\n setting as (-10, 10)\n");
    }

    ALLOC_ARRAY(costs, SIZE, float);
    ALLOC_ARRAY(new_costs, SIZE, float);
    
    if (ga->sel_alg == ROULETTE) {
        ALLOC_ARRAY(probs, SIZE, float);
    } else {
        ALLOC_ARRAY(participants, ga->tour_size, int);
    }
    return 0;

}


static void
free_ga(gaconf *ga)
{
    free(costs);
    free(new_costs);
    free(probs);
    free(participants);
    
    FREE_MATRIX(new_pop, SIZE, 0);
    FREE_MATRIX(pop, SIZE, 1);
    FREE_MATRIX(pairs, SIZE/2, 0);
}


float*
ga(gaconf* ga, float (*func)(float*))
{
    float prev_best = FLT_MAX;
    ga = ga == NULL ? &def : ga;
    int sel_parameter = ga->sel_alg == ROULETTE ? ga->find_max : ga->tour_size;
    int *mm;
    float* ret;

    ASSERT(ga_init(ga) >= 0, NULL,"ga_init fail!\n") 
    srand48(-time(NULL));
    
    generate_initial_pop(ga->ranges);
    calculate_costs(func);

    swap(pop, new_pop, float**);
    swap(costs, new_costs, float*);
    
    for (int it = 0; it < ga->gens; it++) {
        
        sel[ga->sel_alg](sel_parameter);
        crossover_sym();
        mutation(ga->mut_rate);

        calculate_costs(func);
        elitism(ga->elitis * SIZE);
        
        calculate_costs(func); 

        swap(pop, new_pop, float**);
        swap(costs, new_costs, float*);
    }
    mm = max_min(costs);
    puts("");
    swap(pop[0], pop[mm[MAXIMUM^1]], float*);
    return pop[0];
}
