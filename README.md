# genetic-real

Genetic algorithm i wrote for personal use in C.

Written in suckless coding style.
Used on GNU/Linux and OpenBSD.

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
    float min;      /* the minimum if the range */
    float max;      /* the maximum of the range */
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
or simply leave it null to use the defaults.
If ga returns `NULL` there has been an error, otherwise, it returns the best solution as an array.

