#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct simanlconfig
{
    float *ranges; /* Ranges of all particles */
    int dims;      /* Dimensions of the target function. */
    int iters;     /* Number of simulated annealing iterations. */
} simanlconfig;

static void simanl_init(simanlconfig *simanl);
static void find_neighbour(float *C, float *N);
float *simmulated_annealing(float (*E)(float *), simanlconfig *simanl)