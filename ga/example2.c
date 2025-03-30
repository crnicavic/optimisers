#include "genetic.h"


float f(float *x)
{
	return (x[0] - 2) * (x[0] - 2) + x[1] * x[1]; 
}

int main(void)
{
	gaconf conf;
	float ranges[4] = {-10, 10, -10, 10};
	conf.ranges = ranges;	
	conf.dims = 2;
	conf.size = 30;
	conf.tour_size = 0;
	conf.mut_rate = 0.25;
	conf.elitis = 0.2;
	conf.gens = 30;

	float *ret = ga(&conf, f);
	printf("%f, %f\n", ret[0], ret[1]);
	printf("%f\n", f(ret));
	return 0;
}
