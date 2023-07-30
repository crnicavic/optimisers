#include "genetic.h"

#define DIMENSIONS 2

float f2(float* x){
    return (x[0]-2) * (x[0]-2) + x[1] * x[1];
}

extern gaconf def;

int main(){
    if (ga_init(&def) < 0) {
        printf("GA INIT FAILED!\n");
        return 1;
    }
    float* winner = ga(NULL, f2);
    printf("Final: %.10f\n", f2(winner));
    return 0;
}
