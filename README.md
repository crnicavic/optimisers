# genetic-real

Genetic algorithm and Particle Swarm optimization written in C.

```c
ga/genetic.h - ga library header
ga/genetic.c - main ga file
ga/example1.c - neural network training
ga.example2.c - simple 2d function
ga/data/processed.cleveland.data - data the network was trained with
pso/pso.h - pso library header
pso/pso.c - main pso file
```

Keep in mind the neural network is pretty simple, and pretty bad, it
is there mainly for demonstration purposes, and because of my own
curiosity.

The data has been downloaded from https://archive.ics.uci.edu/dataset/45/heart+disease
and i hope they won't sue me.

Usage is simple, include `genetic.h` or `pso.h` and call:
```c
float* ga(gaconf* ga, float (*func)(float*));
or
float *pso(float (*f)(float *), psoconf* pso);
```

# gaconf
```c
typedef struct gaconf{
    float *ranges;  /* the range for each dimension */
    int dims;       /* dimensions of the target function */
    int size;       /* population size */
    int tour_size;  /* the size of the n/2 tournament */
    float mut_rate; /* probability of mutation from 0 to 1 */
    float elitis;   /* percentage of units to be taken into next gen*/ 
    int gens;       /* generation count */
    int find_max;   /* 0 if looking for maximum */
    int sel_alg;    /* selection algorithm */
}gaconf;
```
If the `ga` function returns `NULL` there has been an error, otherwise,
it returns the best solution as an array.
```c
//easier way to create the ranges (works in C)
conf.ranges = {{min_D1, max_D1}, {min_D2, max_D2}, ...
```
if `ranges` are set to `NULL`, the default interval is (-10, 10)

# psoconf
```c
typedef struct psoconf
{
    float *ranges; /* Ranges of all particles */
    int dims;      /* Dimensions of the target function. */
    int size;      /* Size of the swarm. */
    int gens;      /* generation count */
} psoconf;
```
Here ranges work a little different, it only takes 2 elements
and applies to all dimensions, if set to NULL, it uses the
default interval which is (-10, 10).
Also pso doesn't have special functionality for specifying
do you want to find the minimum or maximum, because it isn't
necessary, so just multiply your function by -1 or any hack 
you can think of. If you have a clever or funny solution,
please send it to me!
