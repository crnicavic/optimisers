#include "pso.h"

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
    psoconfig config = {ranges, 1, 15, 1000};
    float *sol = pso(f, &config);
    printf("%f", sol[0]);
    return 0;
}