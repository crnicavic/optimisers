#include "pso.h"


static int SIZE;
static int DIMS;

static float **swarm_pos;
static float **swarm_vel;
static float **particle_best;
static float *global_best;

static float omega = 0.9;
static float phi_p = 2.5;
static float phi_g = 0.5;

static psoconf def = {.ranges = NULL, .dims = 1, .size = 10, .gens=100};

static void
generate_initial_pop(float* ranges)
{
    int i = 0, j;
    for (; i < SIZE; i++) {
        j = 0;
        for (; j < DIMS; j++) {
            swarm_pos[i][j] = RANDOM_FLOAT(ranges[0], ranges[1]);
            swarm_vel[i][j] = 0.0;
            particle_best[i][j] = swarm_pos[i][j];
        }
    }
}

static int 
pso_init(float (*f)(float *), psoconf* pso)
{
    DIMS = pso->dims;
    SIZE = pso->size;
    ALLOC_MATRIX(swarm_pos, SIZE, DIMS, float);
    ALLOC_MATRIX(swarm_vel, SIZE, DIMS, float);
    ALLOC_MATRIX(particle_best, SIZE, DIMS, float);
    ALLOC_ARRAY(global_best, DIMS, float);

    if (pso->ranges == NULL) {
        pso->ranges[0] = -10.0;
        pso->ranges[1] = 10.0;
        printf("No ranges specified. Using default, which is (-10.0, 10.0)");
    }
    printf("%f, %f\n", pso->ranges[0], pso->ranges[1]);

    float *init_best = particle_best[0];
    int i = 1;
    for (; i < pso->size; i++) {
        if (f(particle_best[i]) < f(init_best)) {
            init_best = particle_best[i];
        }
    }
    COPY_ARRAY(init_best, global_best, pso->dims);
}

static void 
find_best(float (*f)(float *))
{
    int i = 0;
    for (; i < SIZE; i++) {
        if (f(swarm_pos[i]) >= f(particle_best[i])) {
            continue;
        }
        COPY_ARRAY(swarm_pos[i], particle_best[i], DIMS);
        if (f(particle_best[i]) < f(global_best)) {
            COPY_ARRAY(particle_best[i], global_best, DIMS);
        }
    }
}

static void
update_swarm(float (*f)(float *), float* ranges)
{
    static float r_p, r_g;
    static float lookahead; /* future pos */
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < DIMS; j++) {
            /* a war crime for control freaks but way cleaner */
            r_p = RANDOM_FLOAT(0, 1) * phi_p;
            r_g = RANDOM_FLOAT(0, 1) * phi_g;
            swarm_vel[i][j] = omega * swarm_vel[i][j] ;
            swarm_vel[i][j] += r_p * (particle_best[i][j] - swarm_pos[i][j]); 
            swarm_vel[i][j] += r_g * (global_best[j] - swarm_pos[i][j]);
            lookahead = swarm_pos[i][j] + swarm_vel[i][j]; 
            if (lookahead < ranges[0] || lookahead > ranges[1]) {
                swarm_pos[i][j] = RANDOM_FLOAT(ranges[0], ranges[1]);
                swarm_vel[i][j] = 0.0;
                continue;
            }

            swarm_pos[i][j] = swarm_pos[i][j] + swarm_vel[i][j];
        }
        find_best(f);
    }
}

float *
pso(float (*f)(float *), psoconf* pso)
{
    ASSERT(pso_init(f, pso) >= 0, NULL, "pso init fail!\n")
    generate_initial_pop(pso->ranges);
    srand(-time(NULL));
    
    for (int iter = 0; iter < pso->gens; iter++) {
        // printf("%f %f\n", *global_best, f(global_best));
        update_swarm(f, pso->ranges);
        phi_p = 2.0 * (float)iter / (float)pso->gens + 2.5;
        phi_g = 2.0 * (float)iter / (float)pso->gens + 2.5;
        omega = 0.5 * (float)iter / (float)pso->gens + 0.9;
    }

    return global_best;
}
