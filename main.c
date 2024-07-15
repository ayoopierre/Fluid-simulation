#include <stdlib.h>
#include <stdbool.h>
#include "init_SDL.h"

int main(int argc, char **argv){
    App *app = init_app();
    if(!app)
        return EXIT_FAILURE;

    while(true){
        //update_physics

        if(!update_app(app)){
            break;
        }
    }

    close_app(app);

    return EXIT_SUCCESS;
}