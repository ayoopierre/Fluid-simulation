#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h> //remove all useless includes afterwards
#include "vector.h"

#define HEIGHT 800
#define WIDTH 600

typedef struct cell{
    float pressure;
    float density;
    vector2 v;
} cell; 

typedef struct Fluid{
    int width;
    int height;
    int num_cells;

    float dt;

    float *pressure;

    float *density;
    float *updated_density;

    vector2 *v;
    vector2 *updated_v;

    float *density_sources;

    float current_max_pressure;
    float current_max_density;
    float diffusion_rate;
} Fluid;
// !!! for performance instead of memcpy swap pointers !!!

// Utils
void init_fluid(Fluid *fluid, int width, int height, float diffusinon_rate, float dt);
int check_allocation(Fluid *fluid);

// Simulation
void handle_sources(Fluid *fluid);
void handle_diffusion(Fluid *fluid);
void handle_advection(Fluid *fluid);
