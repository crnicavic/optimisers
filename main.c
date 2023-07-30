#include "genetic.h"

#define DIMENSIONS 2

float f2(float* x){
    return (x[0]-2) * (x[0]-2) + x[1] * x[1];
}

int main(){
    gaconf conf = {-10, 10, 2, 30, 5, 0.3, 2, 0.01, 30, TRUE, ROULETTE};
    float* winner = ga(&conf, f2);
    printf("Final: %.10f\n", f2(winner));
    return 0;
}
