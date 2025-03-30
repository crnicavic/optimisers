#include "genetic.h"
#include <math.h>

float griewank(float *x)
{
	float sum = 0;
	float product = 1;
	for(int i = 0; i < 2; i++) {
		sum += x[i] * x[i];
		product *= cos(x[i] / sqrt(i+1));
	}
	sum = sum / 4000;
	return 1 + sum - product;
}

int main(void)
{
	gaconf conf;
	float ranges[4] = {-100, 100, -100, 100};
	conf.ranges = ranges;	
	conf.dims = 2;
	conf.size = 100;
	conf.mut_rate = 0.50;
	conf.elitis = 0.1;
	conf.gens = 3000;

	float *ret = ga(&conf, griewank);
	printf("%f, %f\n", ret[0], ret[1]);
	printf("%f\n", griewank(ret));
	return 0;
}
