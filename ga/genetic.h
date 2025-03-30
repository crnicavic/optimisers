#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <string.h>

#define MALLOC_FAIL(str) "%s allocation failed!\n", #str

#define ASSERT(_e, f, ...)  if (!(_e)) { \
                                fprintf(stderr, __VA_ARGS__); \
                                return f; \
                            }
#define swap(x, y, T) do {T swap = x; x = y; y = swap;} while(0);

#define len(arr) sizeof(arr)/sizeof(arr[0])

#define _POSIX_C_SOURCE 200809L


typedef enum selection
{
    ROULETTE = 0,
    TOURNAMENT
}selection;


//TODO: all of these should be function pointers
typedef struct gaconf{
    float *ranges; 
    int dims;       /* dimensions of the target function */
    int size;       /* population size */
    float mut_rate; /* probability of mutation from 0 to 1 */
    float elitis;   /* percentage of units to be taken into next gen*/ 
    int gens;       /* generation count */
}gaconf;

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
