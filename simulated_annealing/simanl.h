#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct simanlconf
{
    float *ranges; /* Ranges of all particles */
    int dims;      /* Dimensions of the target function. */
    int iters;     /* Number of simulated annealing iterations. */
} simanlconf;

static void simanl_init(simanlconf *simanlc);
static void find_neighbour(float *C, float *N);
float *simmulated_annealing(float (*E)(float *), simanlconf *simanlc)