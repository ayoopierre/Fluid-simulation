#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h> //remove all useless includes afterwards
#include "vector.h"

typedef struct Fluid{
    int height;
    int width;
    int num_cells;

    float dt;
    float diffusion_rate;
    float viscosity;

    float *v_x;
    float *v_y;

    float *updated_v_x;
    float *updated_v_y;

    float *density;
    float *updated_density;
} Fluid;

//Utils
void init_fluid(Fluid *fluid, int height, int width, float diff, float visc, float dt);
int check_alloc(Fluid *fluid);
void swap_pointers(void **a, void **b);
void lin_solve(int b, float *updated_x, float *x, float a, float c, int iter, int width, int height);

//Simulation
void add_density(Fluid *fluid, int x, int y, float amount);
void add_velocity(Fluid *fluid, int x, int y, float a_x, float a_y);
void handle_diffusion(int b, int width, int height, float *updated_x, float *x, float diff, float dt, int iterations);
void handle_projection(int width, int height, float *v_x, float *v_y, float *p, float *div, int iterations);
void handle_advection(int b, int width, int height, float *updated_d, float *d, float *v_X, float *v_y, float dt);
void set_bnd(int b, float *x, int width, int height);

void update_fluid_state(Fluid *fluid);