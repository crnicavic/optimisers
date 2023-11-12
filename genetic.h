#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#define MALLOC_FAIL(str) "%s allocation failed!\n", #str

#define ASSERT(_e, f, ...)  if (!(_e)) { \
                                fprintf(stderr, __VA_ARGS__); \
                                return f; \
                            }
/* v- variable(name) s-size t-type */
#define ALLOC_ARRAY(v, s, t) v = (t*) malloc(sizeof(t) * s); \
                        ASSERT(v != NULL, FAIL, MALLOC_FAIL(v))

/* v- variable(name) r-rows c-cols t-type */
#define ALLOC_MATRIX(v, r, c, t) do { int rr = 0;    \
                            v = (t**) malloc(sizeof(t*) * r); \
                            ASSERT(v != NULL, FAIL, MALLOC_FAIL(v)); \
                            for(;rr < r; ++rr) { \
                            v[rr] = (t*) malloc(sizeof(t) * c); \
                            ASSERT(v[rr] != NULL, FAIL, MALLOC_FAIL(v)); } \
                            } while(0);

/* v-variable r-rows s-row to start from */
#define FREE_MATRIX(v, r, s) do { int rr = s; \
                                for(; rr < r; rr++) { \
                                    free(v[rr]); \
                                } \
                                if(!s) { free(v); } } while(0);

#define swap(x, y, T) do {T swap = x; x = y; y = swap;} while(0);

#define _POSIX_C_SOURCE 200809L
#define TRUE 1
#define FALSE 0
#define FAIL -1


typedef enum selection
{
    ROULETTE = 0,
    TOURNAMENT
}selection;

typedef struct gaconf{
    float *ranges; 
    int dims;       /* dimensions of the target function */
    int size;       /* population size */
    int tour_size;  /* the size of the n/2 tournament */
    float mut_rate; /* probability of mutation from 0 to 1 */
    float elitis;   /* percentage of units to be taken into next gen*/ 
    int gens;       /* generation count */
    int find_max;   /* 0 if looking for maximum */
    int sel_alg;    /* selection algorithm */
}gaconf;

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
