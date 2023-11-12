#include "genetic.h"

#define DIMENSIONS 2

float f2(float* x){
        return (x[0]-2) * (x[0]-2) + x[1] * x[1];
}

extern gaconf def;

int main(){
        float test[2][2] = {{-1.0, 3.0}, {-2.0, 3.0}};
        gaconf conf = {test, 2, 30, 5, 0.3, 0.1, 30, FALSE, ROULETTE};
        float* winner = ga(&conf, f2);
        printf("Final: %.10f\n", f2(winner));
        return 0;
}
