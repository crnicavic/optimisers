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
static void crossover_sym();
static void mutation(float mut_rate);
static void elitism(int elitis);
static void** ndarr(const int* shape, int dims, size_t element_size);
static int ga_init(gaconf *ga);
static void free_ga(gaconf *ga);
float* ga(gaconf* ga, float (*func)(float*));

static int DIMS;
static int SIZE;

/* an array that contains all the indicies */
static int *idx;
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
	int it = 0;
	float sum_probs = 0;

	for(; it < SIZE; it++) sum_probs += costs[it];
	
	it = 0;
	for(; it < SIZE; it++) probs[it] = costs[it] / sum_probs;

	it = 0;
    for (; pair < SIZE/2; pair++) {
        crd = 0;
        for (; crd < 2; crd++) {
            pairs[pair][crd] = spin();
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
            new_pop[chr][crd] += r * pop[pairs[pair][1]][crd];
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


/* qsort callback */
static int
compare_costs(const void* a, const void* b)
{
	return costs[*(int*)a] - costs[*(int*)b];
}


static inline void 
elitism(int elitis)
{
	qsort(idx, SIZE, sizeof(int), compare_costs);
}


static void** ndarr
(const int* shape, int dims, size_t element_size)
{
	int i = 0;
	int k = 0;
	size_t total_elements = shape[0];
	size_t total_pointers = shape[dims-2];
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

	level_size = total_elements / shape[dims-1];
	k = total_pointers - level_size;
	/* set in-data pointers */
	while(k < total_pointers) {
		base_ptr[k] = data_iterator;
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
	
	idx = (int*) ndarr(&ga->size, 1, sizeof(int));
	for (int i = 0; i < SIZE; i++) idx[i] = i;

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
    float* ret = malloc(ga->size * sizeof(float));
	int n = 20;

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
        //elitism(ga->elitis * SIZE);
        
        swap(pop, new_pop, float**);
        swap(costs, new_costs, float*);
    }
	memcpy(ret, pop[idx[0]], ga->size * sizeof(float));
	free_ga(ga);
    return ret;
}
