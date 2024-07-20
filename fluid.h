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

    float *v_x;
    float *updated_v_x;

    float *v_y;
    float *updated_v_y;

    float *density_sources;

    float current_max_pressure;
    float current_max_density;
    float diffusion_rate;
    float viscosity;
} Fluid;


// Utils
void init_fluid(Fluid *fluid, int width, int height, float diffusinon_rate, float viscosity, float dt);
void delete_fluid(Fluid *fluid);
int check_allocation(Fluid *fluid);
void randomize_density(Fluid *fluid, float max_density);
void randomize_velocity_field(Fluid *fluid, float min_v, float max_v);
void swap_pointers(float **a, float **b);

// Simulation
void handle_sources(Fluid *fluid);
void handle_diffusion(int width, int height, float *arr, float *updated_arr, float diffusion_rate, float dt);
void handle_advection(int width, int height, float *arr, float *updated_arr, float *v_x, float *v_y, float dt);
void handle_incompresibility(int width, int height, float *v_x, float *v_y, float *updated_v_x, float *updated_v_y);

void update_density(Fluid *fluid);
void update_velocity_field(Fluid *fluid);

void update_fluid_state(Fluid *fluid);

//current to do add velocity field and updating max_density