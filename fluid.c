#include "fluid.h"

void init_fluid(Fluid *fluid, int height, int width, float diff, float visc, float dt){
    fluid->height = height;
    fluid->width = width;
    fluid->num_cells = width * height;

    fluid->diffusion_rate = diff;
    fluid->viscosity = visc ;
    fluid->dt = dt;

    fluid->v_x = calloc(fluid->num_cells, sizeof(float));
    fluid->v_y = calloc(fluid->num_cells, sizeof(float));
    fluid->updated_v_x = calloc(fluid->num_cells, sizeof(float));
    fluid->updated_v_y = calloc(fluid->num_cells, sizeof(float));
    fluid->density = calloc(fluid->num_cells, sizeof(float));
    fluid->updated_density = calloc(fluid->num_cells, sizeof(float));
}

int check_alloc(Fluid *fluid){
    if(NULL == fluid->v_x) return 1;
    if(NULL == fluid->v_y) return 1;
    if(NULL == fluid->updated_v_x) return 1;
    if(NULL == fluid->updated_v_y) return 1;
    if(NULL == fluid->density) return 1;
    if(NULL == fluid->updated_density) return 1;
    return 0;
}

void swap_pointers(void **a, void **b){
    void *temp = *a;
    *a = *b;
    *b = temp;
}

void lin_solve(int b, float *updated_x, float *x, float a, float c, int iter, int width, int height){
    // Implemenatation of Gauss-Seidel method
    float inv_c = 1.0 / c;
    for(int k = 0; k < iter; k++){
        for(int i = 1; i < width - 1; i++){
            for(int j = 1; j < height - 1; j++){
                updated_x[i + width * j] = (x[i + j * width] + a * (
                    updated_x[i - 1 + j * width] +
                    updated_x[i + 1 + j * width] +
                    updated_x[i + (j - 1) * width] + 
                    updated_x[i + (j + 1) * width]
                )) * inv_c;
            }
        }
        // handle bounds;
        set_bnd(b, x, width, height);
    }
}

// Simulaiton

void add_density(Fluid *fluid, int x, int y, float amount){
    fluid->density[x + fluid->width * y] += amount;
    fluid->updated_density[x + fluid->width * y] += amount;
}

void add_velocity(Fluid *fluid, int x, int y, float a_x, float a_y){
    fluid->v_x[x + fluid->width * y] += a_x;
    fluid->v_y[x + fluid->width * y] += a_y;
    fluid->updated_v_x[x + fluid->width * y] += a_x;
    fluid->updated_v_y[x + fluid->width * y] += a_y;
}

void handle_diffusion(int b, int width, int height, float *updated_x, float *x, float diff, float dt, int iterations){
    float a = dt * diff * (width - 2) * (height - 2);
    lin_solve(b, updated_x, x, a, 1 + 6 * a, iterations, width, height);
}

void handle_projection(int width, int height, float *v_x, float *v_y, float *p, float *div, int iterations){
    for(int i = 1; i < width - 1; i++){
        for(int j = 1; j < height - 1; j++){
            div[i + j * width] = -0.5f * (
                    v_x[i + 1 + j * width]
                  - v_x[i - 1 + j * width]
                  + v_y[i + (j + 1) * width] 
                  - v_y[i + (j - 1) * width]
            ) / height;
            p[i + j * width] = 0;
        }
    }
    // handle bounds for p and div
    set_bnd(0, div, width, height);
    set_bnd(0, p, width, height);
    lin_solve(0, p, div, 1, 6, iterations, width, height);

    for(int i = 1; i < width - 1; i++){
        for(int j = 1; j < height - 1; j++){
            v_x[i + j * width] -= 0.5f * (p[i + 1 + j * width] - p[i - 1 + j * width]) * width;
            v_y[i + j * width] -= 0.5f * (p[i + (j + 1) * width] - p[i + (j - 1) * width]) * height;
        }
    }
    //handle bounds for v_x  and v_y
    set_bnd(1, v_x, width, height);
    set_bnd(2, v_y, width, height);
}

void handle_advection(int b, int width, int height, float *updated_d, float *d, float *v_x, float *v_y, float dt){
    float i0, i1, j0, j1;
    float dtx = dt * (width - 2);
    float dty = dt * (height - 2);
    float tmp1, tmp2, x, y;
    float s0, s1, t0, t1;
    float ifloat, jfloat;
    float w_f = width;
    float h_f = height;

    for(int i = 1, ifloat = 1; i < width - 1; i++, ifloat++){
        for(int j = 1, jfloat = 1; j < height -1; j++, jfloat++){
            tmp1 = dtx * v_x[i + j * width];
            tmp2 = dty * v_y[i + j * width];

            x    = ifloat - tmp1; 
            y    = jfloat - tmp2;
            
            if(x < 0.5f) x = 0.5f; 
            if(x > w_f + 0.5f) x = w_f + 0.5f; 
            i0 = floorf(x);                 
            i1 = i0 + 1.0f;
            if(y < 0.5f) y = 0.5f; 
            if(y > h_f + 0.5f) y = h_f + 0.5f; 
            j0 = floorf(y);
            j1 = j0 + 1.0f; 
                
            s1 = x - i0; 
            s0 = 1.0f - s1; 
            t1 = y - j0; 
            t0 = 1.0f - t1;
            
            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;
                
            updated_d[i + j * width] =     
                s0 * ( t0 * d[i0i, j0i * width]
                    + t1 * d[i0i + j1i * width])
               +s1 * ( t0 * d[i1i + j0i * width]
                    + t1 * d[i1i + j1i * width]);
        }
    }
    // handle
    set_bnd(b, d, width, height);
}

void set_bnd(int b, float *x, int width, int height){
    for(int i = 1; i < width; i++) {
        x[i] = b == 2 ? -x[i + width] : x[i + width];
        x[i + (width - 1) * width] = b == 2 ? -x[i + (width - 2) * width] : x[i + (width - 2) * width];
    }

    for(int j = 1; j < height - 1; j++) {
        x[j * width] = b == 1 ? -x[1 + j * width] : x[1 + j * width];
        x[width - 1 + j * width] = b == 1 ? -x[width - 2 + j * width] : x[width - 2 + j * width];
    }
    
    x[0] = (x[1] + x[width]) / 2.0f;
    x[width - 1] = (x[width - 2] + x[2 * width - 1]) / 2.0f;
    x[width * (height - 1)] = (x[width - 2] + x[2 * width - 1]) / 2.0f;
    x[width * height - 1] = (x[width - 1 + (height - 1) * width] + x[width * height - 2]) / 2.0f;
}

void update_fluid_state(Fluid *fluid){
    handle_diffusion(1, fluid->width, fluid->height, fluid->v_x, fluid->updated_v_x, fluid->viscosity, fluid->dt, 4);
    handle_diffusion(2, fluid->width, fluid->height, fluid->v_y, fluid->updated_v_y, fluid->viscosity, fluid->dt, 4);

    handle_projection(fluid->width, fluid->height, fluid->v_x, fluid->v_y, fluid->updated_v_x, fluid->updated_v_y, 4);

    handle_advection(1, fluid->width, fluid->height, fluid->updated_v_x, fluid->v_x, fluid->v_x, fluid->v_y, fluid->dt);
    handle_advection(2, fluid->width, fluid->height, fluid->updated_v_y, fluid->v_y, fluid->v_x, fluid->v_y, fluid->dt);

    handle_projection(fluid->width, fluid->height, fluid->updated_v_x, fluid->updated_v_y, fluid->v_x, fluid->v_y, 4);

    handle_diffusion(0, fluid->width, fluid->height, fluid->updated_density, fluid->density, fluid->diffusion_rate, fluid->dt, 4);
    handle_advection(0, fluid->width, fluid->height, fluid->density, fluid->updated_density, fluid->updated_v_x, fluid->updated_v_y, fluid->dt);
}