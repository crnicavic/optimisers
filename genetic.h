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
static void generate_initial_pop(float *ranges);
static void calculate_costs(float(*f)(float*));
static int greater(float a, float b);
static int lesser(float a, float b);
int quicksort_partition(int start, int stop);
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

