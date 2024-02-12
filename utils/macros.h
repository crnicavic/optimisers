#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* warning - Warning string */
#define WARN(warning)  \
    puts("WARNING\n"); \
    puts(warning);

#define RANDOM_FLOAT(lower, upper) \
    ((float)rand() / (float)RAND_MAX) * (upper - lower) + lower

/* array1 - Copy from array; array2 - Copy to array; size - Size of arrays (equal) */
#define COPY_ARRAY(array1, array2, size) \
    for (int i = 0; i < size; i++)       \
    {                                    \
        array2[i] = array1[i];           \
    }

#define MALLOC_FAIL(str) "%s allocation failed!\n", #str

#define ASSERT(_e, f, ...)            \
    if (!(_e))                        \
    {                                 \
        fprintf(stderr, __VA_ARGS__); \
    }

/* v- variable(name) s-size t-type */
#define ALLOC_ARRAY(v, s, t)        \
    v = (t *)malloc(sizeof(t) * s); \
    ASSERT(v != NULL, FAIL, MALLOC_FAIL(v))

/* v- variable(name) r-rows c-cols t-type */
#define ALLOC_MATRIX(v, r, c, t)                         \
    do                                                   \
    {                                                    \
        int rr = 0;                                      \
        v = (t **)malloc(sizeof(t *) * r);               \
        ASSERT(v != NULL, FAIL, MALLOC_FAIL(v));         \
        for (; rr < r; ++rr)                             \
        {                                                \
            v[rr] = (t *)malloc(sizeof(t) * c);          \
            ASSERT(v[rr] != NULL, FAIL, MALLOC_FAIL(v)); \
        }                                                \
    } while (0);

/* v-variable r-rows s-row to start from */
#define FREE_MATRIX(v, r, s) \
    do                       \
    {                        \
        int rr = s;          \
        for (; rr < r; rr++) \
        {                    \
            free(v[rr]);     \
        }                    \
        if (!s)              \
        {                    \
            free(v);         \
        }                    \
    } while (0);

#define swap(x, y, T) \
    do                \
    {                 \
        T swap = x;   \
        x = y;        \
        y = swap;     \
    } while (0);

#define len(arr) sizeof(arr) / sizeof(arr[0])

#define _POSIX_C_SOURCE 200809L
#define TRUE 1
#define FALSE 0
#define FAIL -1