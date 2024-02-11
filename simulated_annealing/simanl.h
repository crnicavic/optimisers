#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct simanlconf
{
    float *ranges; /* Ranges of all particles */
    int dims;      /* Dimensions of the target function. */
    int iters;     /* Number of simulated annealing iterations. */
} simanlconf;

static void simanl_init(simanlconf *simanl);
static void find_neighbour(float *C, float *N, simanlconf *simanl);
float *simmulated_annealing(float (*E)(float *), simanlconf *simanl);