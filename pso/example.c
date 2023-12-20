#include "pso.h"
#include <math.h>

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
    float ranges[] = {1.5, 6.0};
    float *sol = pso(f, 1, 15, ranges, 1000);
    printf("%f", sol[0]);
    return 0;
}