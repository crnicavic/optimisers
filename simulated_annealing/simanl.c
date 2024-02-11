#include "simanl.h"
#include "../utils/macros.h"

static int DIMS;
static float *C = NULL, *N = NULL;
static simanlconfig *config;

static void simanl_init(simanlconfig *simanl)
{
    srand(-time(NULL));
    config = simanl;

    if (config->ranges == NULL)
    {
        ALLOC_ARRAY(simanl->ranges, 2, float);
        config->ranges[0] = -10.0;
        config->ranges[1] = 10.0;
        WARN("No ranges specified. Using default, which is (-10.0, 10.0)\n");
    }

    DIMS = config->dims;

    ALLOC_ARRAY(C, DIMS, float);
    ALLOC_ARRAY(N, DIMS, float);
    for (int i = 0; i < DIMS; i++)
    {
        C[i] = RANDOM_FLOAT(config->ranges[0], config->ranges[1]);
    }
}

static void find_neighbour(float *C, float *N)
{
    for (int i = 0; i < DIMS; i++)
    {
        N[i] = C[i] + RANDOM_FLOAT(config->ranges[0] / 2.0, config->ranges[1] / 2.0); // Should be discussed later on.
        if (N[i] < config->ranges[0] || N[i] > config->ranges[1])
        {
            /* If out of space range, generate a random solution. */
            N[i] = RANDOM_FLOAT(config->ranges[0], config->ranges[1]);
        }
    }
}

float *simmulated_annealing(float (*E)(float *), simanlconfig *simanl)
{
    simanl_init(simanl);
    int ITERS = simanl->iters;

    /* Minimization problem */
    for (int iter = 0; iter < ITERS; iter++)
    {
        float T = 1.0 - ((float)iter + 1.0) / (float)ITERS;
        find_neighbour(C, N);
        float deltaE = E(N) - E(C);
        if (deltaE < 0 || T > RANDOM_FLOAT(0, 1))
        {
            COPY_ARRAY(N, C, DIMS);
        }
    }
    free(N);
    return C;
}