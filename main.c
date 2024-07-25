#include <stdlib.h>
#include <stdbool.h>
#include "init_SDL.h"

int main(int argc, char **argv){
    App *app = init_app();
    if(!app) return EXIT_FAILURE;

    Fluid fluid;
    init_fluid(&fluid, 60, 80, 0.05f, 3.0f, 0.001);
    
    if(check_alloc(&fluid)) return EXIT_FAILURE;

    while(true){
        //update_physics
        //update_fluid_state(&fluid);
        if(!update_app(app, &fluid)){
            break;
        }
    }

    close_app(app);

    return EXIT_SUCCESS;
}