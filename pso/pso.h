#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define WARNING(str) \
    puts("WARNING"); \
    puts(str);

#define ALLOC_ARRAY(var, size, type) \
    var = (type *)malloc(size * sizeof(type));

#define ALLOC_MATRIX(var, rows, cols, type)       \
    var = (type **)malloc(rows * sizeof(type *)); \
    for (int i = 0; i < rows; i++)                \
    {                                             \
        ALLOC_ARRAY(var[i], cols, type);          \
    }

#define RANDOM_FLOAT(lower, upper) \
    ((float)rand() / (float)RAND_MAX) * (upper - lower) + lower;

#define COPY_ARRAY(array1, array2, size) \
    for (int i = 0; i < size; i++)       \
    {                                    \
        array2[i] = array1[i];           \
    }

typedef struct psoconf
{
    float *ranges; // Ranges of all particles
    int dims;      // Dimensions of the target function.
    int size;      // Size of the swarm.

} psoconf;

float *pso(float (*f)(float *), int dims, int swarm_size, float *ranges, int iterations);
