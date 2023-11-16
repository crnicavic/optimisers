# genetic-real

Genetic algorithm ~~I wrote for personal use~~ in C.

```c
genetic.h - the main file
example1.c - neural network training
example2.c - simple 2d function
data/processed.cleveland.data - data the network was trained with
```

Keep in mind the neural network is pretty simple, and pretty bad, it
is there mainly for demonstration purposes, and because of my own
curiosity.

The data has been downloaded from https://archive.ics.uci.edu/dataset/45/heart+disease
and i hope they won't sue me.

Usage is simple, include `genetic.h` and call:
```c
static float* ga(gaconf* ga, float (*func)(float*));
```
where gaconf is defined as:
```c
typedef struct gaconf{
    float *ranges; 
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
