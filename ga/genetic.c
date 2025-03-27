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
static int partition(int start, int stop);
static int find_kth(int start, int stop, int k);
static int* max_min(float *arr);
static void prob_min(int* max_min);
static inline int spin();
static void roulette();
static int tournament(int participant_count); 
static void brackets(int participant_count);
static void crossover_sym();
static void mutation(float mut_rate);
static void elitism(int elitis);
static void** ndarr(int* shape, int dims, size_t element_size);
static int ga_init(gaconf *ga);
static void free_ga(gaconf *ga);
float* ga(gaconf* ga, float (*func)(float*));

static int DIMS;
static int SIZE;
static int MAXIMUM;

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
    int chr = 0, crd = 0;
    float range, min; 

    for (; chr < SIZE; chr++){
        crd = 0;
        range = ranges[chr * 2 + 1] - ranges[chr * 2];
        for (; crd < DIMS; crd++){
            new_pop[chr][crd] = (drand48() * range) + ranges[chr*2];
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
partition(int start, int stop)
{
    int pivot = start;
    int it = start;
    for (; it < stop; it++) {
        if (new_costs[stop] > new_costs[it]) {
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

/* finds both the max and the min */
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


static void
prob_min(int* max_min)
{
    int chr=0;
    float range = costs[max_min[0]] - costs[max_min[1]];
    float min = costs[max_min[1]];
    for (; chr < SIZE; chr++) {
        probs[chr] = 1 - (costs[chr] - min) / range;    
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
roulette()
{
    int chr = 0, crd = 0;
    int pair = 0;
    int *r = max_min(costs); 

    prob_min(r);

    for (; pair < SIZE/2; pair++) {
        crd = 0;
        for (; crd < 2; crd++) {
            pairs[pair][crd] = spin();
        }
    }
}


static void
brackets(int participant_count)
{
    int parent = 0;
    int pair = 0;
    for (; pair < SIZE/2; pair++) {
        parent = 0;
        for (; parent < 2; parent++) {
            pairs[pair][parent] = tournament(participant_count);
        }
    }
}


/* symmetrical crossover */
static void
crossover_sym()
{
    int chr = 0, crd = 0, pair = 0;
    float r;

    while (chr < SIZE) {
        r = (float) drand48();

        for (; crd < DIMS; crd++) {
            new_pop[chr][crd] = r * pop[pairs[pair][0]][crd];
            new_pop[chr][crd] += (1-r) * pop[pairs[pair][1]][crd];
        }
        
        crd = 0;
        chr++;

        for (; crd < DIMS; crd++) {
            new_pop[chr][crd] = (1-r) * pop[pairs[pair][0]][crd]; 
            new_pop[chr][crd] =+ r * pop[pairs[pair][1]][crd];
        }

        chr++;
        pair++;
    }
}


static void
mutation(float mut_rate)
{
    float r;
    int chr = 0, crd = 0;
    for (; chr < SIZE; chr++) {
        crd = 0;
        for (; crd < DIMS; crd++) {  
            r = drand48() < mut_rate ? (drand48() * 2) - 1 : 0;
            new_pop[chr][crd] += r;
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


void** ndarr
(int* shape, int dims, size_t element_size)
{
	int i = 0;
	int k = 0;
	size_t total_elements = shape[0];
	size_t total_pointers = shape[dims-2];
	size_t total_length = 0;
	size_t total_size = 0;
	size_t level_size = 0;

	for(i = 1; i < dims; i++) {
		total_elements *= shape[i];
	}
	
	for(i = dims-3; i >= 0; i--) {
		total_pointers *= shape[i];
		total_pointers += shape[i];
	}
	
	total_size = total_elements * element_size + total_pointers * sizeof(void**);
	void** base_ptr = (void**) malloc(total_size);
	if (base_ptr == NULL)
	{
		printf("memory alloc failed!");
		exit(1);
	}

	void **base_data = base_ptr + total_pointers;
	
	void *data_iterator = base_data;
	i = total_pointers - 1;

	level_size = total_elements / shape[dims-1];
	k = total_pointers - level_size;
	int loop_counter = 0; /* debugging purposes, will be removed */
	/* set in-data pointers */
	while(k < total_pointers) {
		base_ptr[k] = data_iterator;
		loop_counter++;
		k++;
		data_iterator += element_size * shape[dims-1];
	}

	int start;
	int end = total_pointers;
	void *pointer_iterator;
	for(i = dims - 3; i >= 0; --i) {
		end -= level_size;
		level_size = level_size / shape[i+1];
		start = end - level_size; 
		pointer_iterator = base_ptr + end;
		k = start;
		printf("start: %d\n", start);
		for(k = start; k < end; k++)
		{
			base_ptr[k] = pointer_iterator;
			pointer_iterator += shape[i+1] * sizeof(void*);
		}
	}
	return base_ptr;
}


/* dumb function that calls a bunch of mallocs */
static int
ga_init(gaconf *ga)
{
    int pair = 0, it = 0;
    float min, max;
 
    SIZE = ga->size;
    DIMS = ga->dims;
	int pop_shape[2] = {SIZE, DIMS};
	int pairs_shape[2] = {SIZE/2, 2};
	

	pop = (float**) ndarr(pop_shape, 2, sizeof(float)); 
	new_pop = (float**) ndarr(pop_shape, 2, sizeof(float));
	pairs = (int**) ndarr(pairs_shape, 2, sizeof(int));
    
    if (ga->ranges == NULL || len(ga->ranges) != DIMS * 2) {
        min = -10;
        max = 10;
        if (ga->ranges != NULL && len(ga->ranges) == 2) {
            min = ga->ranges[0];
            max = ga->ranges[1];
        }

		int ranges_shape = DIMS * 2;
		ga->ranges = (float*) ndarr(&ranges_shape, 1, sizeof(float));
        it = 0;
        for (; it < DIMS; it++){
            ga->ranges[it * 2] = min;
            ga->ranges[it * 2 + 1] = max;
        }
    }

	costs = (float*) ndarr(&ga->size, 1, sizeof(float));
	new_costs = (float*) ndarr(&ga->size, 1, sizeof(float));

	probs = (float*) ndarr(&ga->size, 1, sizeof(float));
    return 0;

}

static void
free_ga(gaconf *ga)
{
    free(costs);
    free(new_costs);
    free(probs);
    free(participants);
    
	free(pop);
	free(new_pop);
	free(pairs);
}

float*
ga(gaconf* ga, float (*func)(float*))
{
    float prev_best = FLT_MAX;
    int *mm;
    float* ret = malloc(ga->size * sizeof(float));
	int n = 20;
	int arr[n];

    ASSERT(ga_init(ga) >= 0, NULL,"ga_init fail!\n") 
    srand48(-time(NULL));
    
    generate_initial_pop(ga->ranges);
    calculate_costs(func);

    swap(pop, new_pop, float**);
    swap(costs, new_costs, float*);
    
    for (int it = 0; it < ga->gens; it++) {
        
		roulette();
        crossover_sym();
        mutation(ga->mut_rate);

        calculate_costs(func);
        elitism(ga->elitis * SIZE);
        
        calculate_costs(func); 

        swap(pop, new_pop, float**);
        swap(costs, new_costs, float*);
    }
    mm = max_min(costs);
    swap(pop[0], pop[mm[1]], float*);
	memcpy(ret, pop[0], ga->size * sizeof(float));
	free_ga(ga);
    return ret;
}
