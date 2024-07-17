# include "fluid.h"

// Utilities
void init_fluid(Fluid *fluid, int width, int height, float diffusinon_rate, float viscosity, float dt){
    fluid->height = height;
    fluid->width = width;
    fluid->diffusion_rate = diffusinon_rate;
    fluid->viscosity = viscosity;
    fluid->dt = dt;

    //fluid->num_cells = (fluid->height + 2) * (fluid->width + 2);
    fluid->num_cells = fluid->height * fluid->width;

    fluid->pressure = calloc(fluid->num_cells, sizeof(float));
    fluid->density = calloc(fluid->num_cells, sizeof(float));
    fluid->updated_density = calloc(fluid->num_cells, sizeof(float));
    fluid->density_sources = calloc(fluid->num_cells, sizeof(float));

    fluid->v_x = calloc(fluid->num_cells, sizeof(float));
    fluid->updated_v_x = calloc(fluid->num_cells, sizeof(float));

    fluid->v_y = calloc(fluid->num_cells, sizeof(float));
    fluid->updated_v_y = calloc(fluid->num_cells, sizeof(float));
}

void delete_fluid(Fluid *fluid){
    free(fluid->density);
    free(fluid->pressure);
    free(fluid->updated_density);
    free(fluid->v_x);
    free(fluid->updated_v_x);
    free(fluid->v_y);
    free(fluid->updated_v_y);
    free(fluid->density_sources);
}

int check_allocation(Fluid *fluid){
    if(NULL == fluid->pressure) return 1;
    if(NULL == fluid->density) return 2;
    if(NULL == fluid->updated_density) return 3;
    if(NULL == fluid->density_sources) return 4;
    if(NULL == fluid->v_x) return 5;
    if(NULL == fluid->updated_v_x) return 6;
    if(NULL == fluid->v_y) return 7;
    if(NULL == fluid->updated_v_y) return 8;
    return 0;
}

void randomize_density(Fluid *fluid, float max_density){
    for(int i = 0; i < fluid->num_cells; i++){
        fluid->density[i] = (float)rand()/(float)(RAND_MAX/max_density);
    }
    fluid->current_max_density = max_density;
}

void randomize_velocity_field(Fluid *fluid, float min_v, float max_v){
    for(int i = 0; i < fluid->num_cells; i++){
        fluid->v_x[i] = (float)rand()/(float)(RAND_MAX/(max_v - min_v)) + min_v;
        fluid->v_y[i] = (float)rand()/(float)(RAND_MAX/(max_v - min_v)) + min_v;
    }
}

void swap_pointers(float **a, float **b){
    float *temp = *a;
    *a = *b;
    *b = temp;
}

// Simulation

void handle_sources(Fluid *fluid){
    for(int i = 0; i < fluid->num_cells; i++){
        fluid->density[i] += fluid->density_sources[i] * fluid->dt;
    }
}


void handle_diffusion(int width, int height, float *arr, float *updated_arr, float diffusion_rate, float dt){
    float a = diffusion_rate * dt * width * height;

    for(int k = 0; k < 20; k++){
        for(int x = 1; x < width-1; x++){
            for(int y = 1; y < height-1; y++){
                updated_arr[x + y * width] = (arr[x + y * width] + a * (arr[x - 1 + y * width] + 
                                            arr[x + 1 + y * width] + arr[x + (y - 1) * width] +
                                             arr[x + (y + 1) * width])) / (1 + 4 * a);
                //if(fluid->current_max_density < fluid->updated_density[x + y * w]) fluid->current_max_density = fluid->updated_density[x + y * w];
            }
        }
        // handle bounds;
    }
}

void handle_advection(int width, int height, float *arr, float *updated_arr, float *v_x, float *v_y, float dt){
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * width;
    
    for(i = 1; i < width-1; i++){
        for(j = 1; j < height-1; j++){
            x = i - dt0 * v_x[i + j * width];
            y = j - dt0 * v_y[i + j * width];

            if(x < 0.5) x = 0.5;
            if(x > width + 0.5) x = width + 0.5;
            i0 = (int)x;
            i1 = i0 + 1;

            if(y < 0.5) y = 0.5;
            if(y > height + 0.5) x = height + 0.5;
            j0 = (int)x;
            j1 = j0 + 1;

            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;
            updated_arr[i + j * width] = s0 * (t0 * arr[i0 + j0 * width] + t1 * arr[i0 + j1 * width] + t0 * arr[i1 + j0 * width] + t1 * arr[i1 + j0 * width]);
            //if(fluid->updated_density[i + j * w] > fluid->current_max_density) fluid->current_max_density = fluid->updated_density[i + j * w];
        }
    }
    // handle bounds;
}

void handle_incompresibility(int width, int height, float *v_x, float *v_y, float *updated_v_x, float *updated_v_y){
    for(int i = 1; i < width - 1; i++){
        for(int j = 1; j < height - 1; j++){
            updated_v_y[i + j * width] = -0.5 * 1 / width * (v_x[i + 1 + j * width] - v_x[i - 1 + j * width] + v_y[i + (j + 1) * width] - v_y[i + (j - 1) * width]);
            updated_v_x[i + j * width] = 0;
        }
    }

    for(int k = 0; k < 20; k++){
        for(int i = 1; i < width - 1; i++){
            for(int j = 1; j < height - 1; j++){
                updated_v_x[i + j * width] = (updated_v_y[i - 1 + j *width] + updated_v_x[i + 1 + j * width] + updated_v_x[i + (j + 1) * width] + updated_v_x[i + (j - 1) * width]) / 4;
            }
        }
    }

    for(int i = 1; i < width - 1; i++){
        for(int j = 1; j < height - 1; j++){
            v_x[i + j * width] -= 0.5 * (updated_v_x[i + 1 + j * width] - updated_v_x[i - 1 + j * width]);
            v_y[i + j * width] -= 0.5 * (updated_v_x[i + (j + 1) * width] - updated_v_x[i + (j - 1) * width]);
        }
    }
}


void update_density(Fluid *fluid){
    //sources
    handle_diffusion(fluid->width, fluid->height, fluid->density, fluid->updated_density, fluid->diffusion_rate, fluid->dt);
    swap_pointers(&fluid->density, &fluid->updated_density);
    handle_advection(fluid->width, fluid->height, fluid->density, fluid->updated_density, fluid->v_x, fluid->v_y, fluid->dt);
    swap_pointers(&fluid->density, &fluid->updated_density);
}

void update_velocity_field(Fluid *fluid){
    handle_diffusion(fluid->width, fluid->height, fluid->v_x, fluid->updated_v_x, fluid->viscosity, fluid->dt);
    swap_pointers(&fluid->v_x, &fluid->updated_v_x);
    handle_diffusion(fluid->width, fluid->height, fluid->v_y, fluid->updated_v_y, fluid->viscosity, fluid->dt);
    swap_pointers(&fluid->v_y, &fluid->updated_v_y);
    handle_incompresibility(fluid->width, fluid->height, fluid->v_x, fluid->v_y, fluid->updated_v_x, fluid->updated_v_y);
    swap_pointers(&fluid->v_y, &fluid->updated_v_y);
    swap_pointers(&fluid->v_x, &fluid->updated_v_x);
    handle_advection(fluid->width, fluid->height, fluid->v_x, fluid->updated_v_x, fluid->v_x, fluid->v_y, fluid->dt);
    handle_advection(fluid->width, fluid->height, fluid->v_y, fluid->updated_v_y, fluid->v_x, fluid->v_y, fluid->dt);
    handle_incompresibility(fluid->width, fluid->height, fluid->v_x, fluid->v_y, fluid->updated_v_x, fluid->updated_v_y);
}

void update_fluid_state(Fluid *fluid){
    update_velocity_field(fluid);
    update_density(fluid);
}