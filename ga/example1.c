#include "genetic.h"
#include <math.h>

#define DIMENSIONS 14
#define WEIGHT_COUNT 13*16+16
#define NEURON_COUNT 13+16+1
#define PATIENTS 297
#define LAYER_COUNT 4
#define TRAINING_SIZE 200
#define FILE_NAME "data/processed.cleveland.data" 


/* weights is essentialy a matrix
 * where every row has |prev_layer_neurons| columns
 * so between the first and second layer
 * 16 rows (because the second layer has 16 rows)
 * each one having 13 columns
 */
float weights[WEIGHT_COUNT];
float activations[NEURON_COUNT];
float data[PATIENTS][14];

const int npl[] = {13, 16, 1};
int gen = 0;


/* TODO: Sigmoid */

float
sigmoid(float val) {
    return val / (1 + abs(val));
}

float 
func(float* weights) 
{
    int neuron_it = 0;
    int pl_neuron_it = 0; /* previous layer it */
    float temp;
    int it = 0;
    static int training_it = 0;
    for (; neuron_it < npl[0]; neuron_it++) {
        activations[neuron_it] = sigmoid(data[training_it][neuron_it]);
    }
    neuron_it = npl[0];
    for (; neuron_it <= npl[0] + npl[1]; neuron_it++) {
        activations[neuron_it] = 0;
    }
    neuron_it = npl[0];
    for (; neuron_it < npl[0] + npl[1]; neuron_it++) {
        /* iterate through weights, somehow?? */
        /* do all the sums and stuff */
        pl_neuron_it = 0;
        temp = 0;
        for (; pl_neuron_it < npl[0]; pl_neuron_it++) {
            temp = activations[pl_neuron_it];
            temp *= weights[neuron_it * npl[1] + pl_neuron_it];
            activations[neuron_it] += sigmoid(temp);
        }
        activations[neuron_it] = sigmoid(activations[neuron_it]);
    }
    neuron_it = npl[0] + npl[1];
    pl_neuron_it = npl[0] + npl[1]-2;
    temp = 0;
    for (; pl_neuron_it > npl[1] ; pl_neuron_it--) {
        temp = activations[pl_neuron_it];
        temp *= weights[WEIGHT_COUNT-1-pl_neuron_it+npl[0]+npl[1]]; 
        activations[neuron_it] += sigmoid(temp);
    }
    activations[neuron_it] = sigmoid(activations[neuron_it]) * 4;

    training_it++;
    training_it = training_it >= TRAINING_SIZE ? training_it : 0;
    
    temp = floor(activations[neuron_it]);
    temp = temp < 0 ? 0 : temp;
    temp = temp > 4 ? 4 : temp;
    
    activations[neuron_it] = temp;
/*
    for (it = 0; it < NEURON_COUNT; it++) {
        printf("%.2f, ", activations[it]);
    }
*/
    return activations[neuron_it];
}

float
f2(float *weights)
{
    /* call func() for every piece of training data */
    float err_sum = 0;
    int it = 0;
    for (; it < TRAINING_SIZE; it++) {
        //err_sum += powf(data[it][13] - func(weights), 2);
        err_sum += data[it][13] == func(weights);
    }
    return err_sum;
}

int 
main()
{
    char temp[70];
    int file_it = 0;
    int it = 0;
    int col = 0;
    int row = 0;
    int buf[2];
    int parse_dec = 0;
    float* winner;
    float err = 0;

    gaconf conf = {NULL, WEIGHT_COUNT, 100, 5, 0.3, 0.1, 100, TRUE, TOURNAMENT};
    
    FILE *f = fopen(FILE_NAME, "r+");
    ASSERT(f != NULL, FAIL, "%s failed to open", FILE_NAME);

    while (fscanf(f,"%s", temp) != EOF){
        it = 0;
        buf[0] = 0; buf[1] = 0;
        col = 0;
        while (temp[it] != '\0') {
            if(temp[it] >= '0' && temp[it] <= '9') {
                buf[parse_dec] *= 10;
                buf[parse_dec] += temp[it] - '0';
            } else if (temp[it] == '.') {
                buf[1] = 0;
                parse_dec = 1;
            } else if (temp[it] == ',') {
                data[file_it][col] = buf[0] + buf[1] / 10;
                parse_dec = 0;
                buf[0] = 0; buf[1] = 0;
                col++;
            } else {
                printf("unknown action for %c", temp[it]);
                exit(1);
            }
            
            it++;
        }
        if (temp[it] == '\0') {
            data[file_it][col] = buf[0] + buf[1] / 10;
            parse_dec = 0;
            buf[0] = 0; buf[1] = 0;
        }
        file_it++;
    }
    winner = ga(&conf, f2); 
    err = f2(winner);
    printf("%.4f\n", err);
    return 0;
}


