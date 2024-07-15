# include "fluid.h"

// Utilities
void init_fluid(Fluid *fluid, int width, int height, float diffusinon_rate, float dt){
    fluid->height = height;
    fluid->width = width;
    fluid->diffusion_rate = diffusinon_rate;
    fluid->dt = dt;

    fluid->num_cells = (fluid->height + 2) * (fluid->width + 2);

    fluid->pressure = calloc(fluid->num_cells, sizeof(float));
    fluid->density = calloc(fluid->num_cells, sizeof(float));
    fluid->updated_density = calloc(fluid->num_cells, sizeof(float));
    fluid->density_sources = calloc(fluid->num_cells, sizeof(float));

    fluid->v = calloc(fluid->num_cells, sizeof(vector2));
    fluid->updated_v = calloc(fluid->num_cells, sizeof(vector2));
}

void delete_fluid(Fluid *fluid){
    free(fluid->density);
    free(fluid->pressure);
    free(fluid->updated_density);
    free(fluid->v);
    free(fluid->updated_v);
    free(fluid->density_sources);
}

int check_allocation(Fluid *fluid){
    if(NULL == fluid->pressure) return 1;
    if(NULL == fluid->density) return 2;
    if(NULL == fluid->updated_density) return 3;
    if(NULL == fluid->density_sources) return 4;
    if(NULL == fluid->v) return 5;
    if(NULL == fluid->updated_v) return 6;
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
        fluid->v[i].x = (float)rand()/(float)(RAND_MAX/(max_v - min_v)) + min_v;
        fluid->v[i].y = (float)rand()/(float)(RAND_MAX/(max_v - min_v)) + min_v;
    }
}

void swap_pointers(float *a, float *b){
    float *temp = a;
    a = b;
    b = temp;
}

// Simulation


void handle_diffusion(Fluid *fluid){
    float a = fluid->diffusion_rate * fluid->dt * fluid->num_cells;
    int w = fluid->width + 2;

    for(int k = 0; k < 20; k++){
        for(int x = 1; x <= fluid->width; x++){
            for(int y = 1; y <= fluid->height; y++){
                fluid->updated_density[x + y * w] = (fluid->density[x + y * w] + a * (fluid->density[x - 1 + y * w] + 
                                            fluid->density[x + 1 + y * w] + fluid->density[x + (y - 1) * w] +
                                             fluid->density[x + (y + 1) * w])) / (1 + 4 * a);
                //if(fluid->current_max_density < fluid->updated_density[x + y * w]) fluid->current_max_density = fluid->updated_density[x + y * w];
            }
        }
        // handle bounds;
    }
}

void handle_advection(Fluid *fluid){
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = fluid->dt * fluid->width;
    int w = fluid->width + 2;
    int h = fluid->height + 2;
    
    for(i = 1; i <= fluid->width; i++){
        for(j = 1; j <= fluid->height; j++){
            x = i - dt0 * fluid->v[i + j * w].x;
            y = j - dt0 * fluid->v[i + j * w].y;

            if(x < 0.5) x = 0.5;
            if(x > fluid->width + 0.5) x = fluid->width + 0.5;
            i0 = (int)x;
            i1 = i0 + 1;

            if(y < 0.5) y = 0.5;
            if(y > fluid->height + 0.5) x = fluid->height + 0.5;
            j0 = (int)x;
            j1 = j0 + 1;

            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;
            fluid->updated_density[i + j * w] = s0 * (t0 * fluid->density[i0 + j0 * w] + t1 * fluid->density[i0 + j1 * w] + t0 * fluid->density[i1 + j0 * w] + t1 * fluid->density[i1 + j0 * w]);
            //if(fluid->updated_density[i + j * w] > fluid->current_max_density) fluid->current_max_density = fluid->updated_density[i + j * w];
        }
    }
    // handle bounds;
}

void update_density(Fluid *fluid){
    //sources
    handle_diffusion(fluid);
    swap_pointers(fluid->density, fluid->updated_density);
    handle_advection(fluid);
    swap_pointers(fluid->density, fluid->updated_density);
    printf("%f\n", fluid->density[1234]);
}


void update_fluid_state(Fluid *fluid){
    update_density(fluid);
}