#include <stdlib.h>
#include <stdbool.h>
#include "init_SDL.h"

int main(int argc, char **argv){
    App *app = init_app();
    if(!app) return EXIT_FAILURE;

    Fluid fluid;
    init_fluid(&fluid, 80, 60, 0.5, 0.5);
    if(check_allocation(&fluid) != 0) return EXIT_FAILURE;

    randomize_density(&fluid, 5.0);
    randomize_velocity_field(&fluid, -1.0, 2.0);

    while(true){
        //update_physics

        update_fluid_state(&fluid);

        if(!update_app(app, &fluid)){
            break;
        }
    }

    close_app(app);
    delete_fluid(&fluid);

    return EXIT_SUCCESS;
}