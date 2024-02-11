#include "simanl.h"

float f(float *x)
{
    return (float)pow(*x, 3) - 3.0 * (float)pow(*x, 2) - 4.0;
}

float g(float *x)
{
    return (*x + 2) * (*x + 2);
}

int main()
{
    float ranges[] = {0.0, 6.0};
    simanlconf simanlconf = {ranges, 1, 1000};
    float *sol = simmulated_annealing(f, &simanlconf);
    printf("%f", sol[0]);
    free(sol);
    return 0;
}