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
    psoconf psoconf = {ranges, 1, 10, 100};
    float *sol = pso(f, &psoconf);
    printf("%f", sol[0]);
    return 0;
}
