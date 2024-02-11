#include "../utils/macros.h"
typedef struct psoconfig
{
    float *ranges; /* Ranges of all particles */
    int dims;      /* Dimensions of the target function. */
    int size;      /* Size of the swarm. */
    int iters;     /* Number of PSO iterations. */
} psoconfig;

static void generate_initial_swarm();
static void pso_init(float (*f)(float *), psoconfig *psoc);
static void find_best(float (*f)(float *));
float *pso(float (*f)(float *), psoconfig *psoc);