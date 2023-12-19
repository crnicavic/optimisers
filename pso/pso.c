#include "pso.h"

static float **swarm_pos;
static float **swarm_vel;
static float **particle_best;
static float *global_best;

static float omega = 0.9;
static float phi_p = 2.5;
static float phi_g = 0.5;

static psoconf conf = {.ranges = NULL, .dims = 1, .size = 10};

static void
generate_initial_pop()
{
    ALLOC_MATRIX(swarm_pos, conf.size, conf.dims, float);
    ALLOC_MATRIX(swarm_vel, conf.size, conf.dims, float);
    ALLOC_MATRIX(particle_best, conf.size, conf.dims, float);
    ALLOC_ARRAY(global_best, conf.dims, float);
    for (int i = 0; i < conf.size; i++)
    {
        for (int j = 0; j < conf.dims; j++)
        {
            swarm_pos[i][j] = RANDOM_FLOAT(conf.ranges[0], conf.ranges[1]);
            swarm_vel[i][j] = 0.0;
            particle_best[i][j] = swarm_pos[i][j];
        }
    }
}

static void pso_init(float (*f)(float *), int dims, int swarm_size, float *ranges)
{
    srand((unsigned)time(NULL));
    conf.dims = dims;
    conf.size = swarm_size >= 1 ? swarm_size : 15;
    ALLOC_ARRAY(conf.ranges, 2, float);
    if (ranges == NULL)
    {
        conf.ranges[0] = -10.0;
        conf.ranges[1] = 10.0;
        WARNING("No ranges specified. Using default, which is (-10.0, 10.0)");
    }
    else
    {
        conf.ranges[0] = ranges[0];
        conf.ranges[1] = ranges[1];
    }

    generate_initial_pop();

    float *init_best = particle_best[0];
    for (int i = 1; i < conf.size; i++)
    {
        if (f(particle_best[i]) < f(init_best))
        {
            init_best = particle_best[i];
        }
    }
    COPY_ARRAY(init_best, global_best, conf.dims);
}

static void find_best(float (*f)(float *))
{
    for (int i = 0; i < conf.size; i++)
    {
        if (f(swarm_pos[i]) < f(particle_best[i]))
        {
            COPY_ARRAY(swarm_pos[i], particle_best[i], conf.dims);
            if (f(particle_best[i]) < f(global_best))
            {
                COPY_ARRAY(particle_best[i], global_best, conf.dims);
            }
        }
    }
}

float *pso(float (*f)(float *), int dims, int swarm_size, float *ranges, int iterations)
{
    pso_init(f, dims, swarm_size, ranges);

    for (int iter = 0; iter < iterations; iter++)
    {
        // printf("%f %f\n", *global_best, f(global_best));
        for (int i = 0; i < conf.size; i++)
        {
            for (int j = 0; j < conf.dims; j++)
            {
                float r_p = RANDOM_FLOAT(0, 1);
                float r_g = RANDOM_FLOAT(0, 1);
                swarm_vel[i][j] = omega * swarm_vel[i][j] + phi_p * r_p * (particle_best[i][j] - swarm_pos[i][j]) + phi_g * r_g * (global_best[j] - swarm_pos[i][j]);
                if (swarm_pos[i][j] + swarm_vel[i][j] < conf.ranges[0] || swarm_pos[i][j] + swarm_vel[i][j] > conf.ranges[1])
                {
                    swarm_pos[i][j] = RANDOM_FLOAT(conf.ranges[0], conf.ranges[1]);
                    swarm_vel[i][j] = 0.0;
                    continue;
                }

                swarm_pos[i][j] = swarm_pos[i][j] + swarm_vel[i][j];
            }
            find_best(f);
        }
        phi_p = 2.0 * (float)iter / (float)iterations + 2.5;
        phi_g = 2.0 * (float)iter / (float)iterations + 2.5;
        omega = 0.5 * (float)iter / (float)iterations + 0.9;
    }

    return global_best;
}
