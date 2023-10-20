# genetic-real

Genetic algorithm i wrote for personal use in C.

Written in suckless coding style (not completely), but the software still sucks. A lot.

Compiled with gcc, because otherwise drand doesn't work, and I didn't want to bother.

```c
genetic.h - the main file
main.c and Makefile - if you wish to test it yourself
```

The genetic.h file has a lot of functions, but the only 2 of concern are
```c
static int ga_init(gaconf *ga);
static float* ga(gaconf* ga, float (*func)(float*));
```
Before calling ga, which is the function that does everything, call ga_init to allocate everything.
You can do it yourself if for some reason dynamic allocation isn't an option.
