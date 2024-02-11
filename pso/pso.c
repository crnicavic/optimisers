#include "pso.h"

static float **swarm_pos;
static float **swarm_vel;
static float **particle_best;
static float *global_best;

static float omega = 0.9;
static float phi_p = 2.5;
static float phi_g = 0.5;

static psoconf *config;
static int DIMS;
static int SIZE;

static void generate_initial_swarm()
{
    ALLOC_MATRIX(swarm_pos, SIZE, DIMS, float);
    ALLOC_MATRIX(swarm_vel, SIZE, DIMS, float);
    ALLOC_MATRIX(particle_best, SIZE, DIMS, float);
    ALLOC_ARRAY(global_best, DIMS, float);

    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < DIMS; j++)
        {
            swarm_pos[i][j] = RANDOM_FLOAT(config->ranges[0], config->ranges[1]);
            swarm_vel[i][j] = 0.0;
            particle_best[i][j] = swarm_pos[i][j];
        }
    }
}

static void pso_init(float (*f)(float *), psoconf *psoc)
{
    srand(-time(NULL));

    config = psoc;
    DIMS = config->dims;
    SIZE = config->size;

    if (config->ranges == NULL)
    {
        ALLOC_ARRAY(config->ranges, 2, float);
        config->ranges[0] = -10.0;
        config->ranges[1] = 10.0;
        WARN("No ranges specified. Using default, which is (-10.0, 10.0)\n");
    }

    generate_initial_swarm();

    float *init_best = particle_best[0];
    for (int i = 1; i < SIZE; i++)
    {
        if (f(particle_best[i]) < f(init_best))
        {
            init_best = particle_best[i];
        }
    }
    COPY_ARRAY(init_best, global_best, DIMS);
}

static void find_best(float (*f)(float *))
{
    for (int i = 0; i < SIZE; i++)
    {
        if (f(swarm_pos[i]) < f(particle_best[i]))
        {
            COPY_ARRAY(swarm_pos[i], particle_best[i], DIMS);
            if (f(particle_best[i]) < f(global_best))
            {
                COPY_ARRAY(particle_best[i], global_best, DIMS);
            }
        }
    }
}

float *pso(float (*f)(float *), psoconf *psoc)
{
    pso_init(f, psoc);
    int ITERS = config->iters;

    for (int iter = 0; iter < ITERS; iter++)
    {
        for (int i = 0; i < SIZE; i++)
        {
            for (int j = 0; j < DIMS; j++)
            {
                float r_p = RANDOM_FLOAT(0, 1);
                float r_g = RANDOM_FLOAT(0, 1);
                swarm_vel[i][j] = omega * swarm_vel[i][j];
                swarm_vel[i][j] += phi_p * r_p * (particle_best[i][j] - swarm_pos[i][j]);
                swarm_vel[i][j] += phi_g * r_g * (global_best[j] - swarm_pos[i][j]);
                if (swarm_pos[i][j] + swarm_vel[i][j] < config->ranges[0] || swarm_pos[i][j] + swarm_vel[i][j] > config->ranges[1])
                {
                    swarm_pos[i][j] = RANDOM_FLOAT(config->ranges[0], config->ranges[1]);
                    swarm_vel[i][j] = 0.0;
                    continue;
                }

                swarm_pos[i][j] = swarm_pos[i][j] + swarm_vel[i][j];
            }
            find_best(f);
        }
        phi_p = -2.0 * (float)(iter + 1) / (float)ITERS + 2.5;
        phi_g = 2.0 * (float)(iter + 1) / (float)ITERS + 2.5;
        omega = -0.5 * (float)(iter + 1) / (float)ITERS + 0.9;
    }

    return global_best;
}