# genetic-real

Genetic algorithm i wrote for personal use in C.

```c
genetic.h - the main file
main.c and Makefile - if you wish to test it yourself
```

Usage is simple, include `genetic.h` and call:
```c
static float* ga(gaconf* ga, float (*func)(float*));
```
where gaconf is defined as:
```c
typedef struct gaconf{
    float *ranges;  /* 
    int dims;       /* dimensions of the target function */
    int size;       /* population size */
    int tour_size;  /* the size of the n/2 tournament */
    float mut_rate; /* probability of mutation from 0 to 1 */
    int elitis;     /* count of units to conserve to elitis */
    int gens;       /* generation count */
    int find_max;   /* 0 if looking for maximum */
    int sel_alg;    /* selection algorithm */
}gaconf;
```
or simply leave it `NULL` to use the defaults.
If the `ga` function returns `NULL` there has been an error, otherwise,
it returns the best solution as an array.
```c
//easier way to create the ranges (works in C)
conf.ranges = {{min_D1, max_D1}, {min_D2, max_D2}, ...
```
if `ranges` are set to `NULL`, the default interval is (-10, 10)
