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
} cell; //used in initial approach, changed for cache performance

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


// Utils
void init_fluid(Fluid *fluid, int width, int height, float diffusinon_rate, float dt);
void delete_fluid(Fluid *fluid);
int check_allocation(Fluid *fluid);
void randomize_density(Fluid *fluid, float max_density);
void randomize_velocity_field(Fluid *fluid, float min_v, float max_v);
void swap_pointers(float *a, float *b);

// Simulation
void handle_sources(Fluid *fluid);
void handle_diffusion(Fluid *fluid);
void handle_advection(Fluid *fluid);

void update_density(Fluid *fluid);

void update_fluid_state(Fluid *fluid);

//current to do add velocity field and updating max_density