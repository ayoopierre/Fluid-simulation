#include "init_SDL.h"

App* init_app(){
    App *app = malloc(sizeof(App));

    app->fps = FPS;

    if(SDL_Init(SDL_INIT_EVERYTHING) < 0){
        printf("Failed to initialise SDL\n");
        free(app);
        return NULL;
    }

    app->window = SDL_CreateWindow("Window", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
    WINDOW_WIDTH, WINDOW_HEIGHT, 0);

    if(!app->window){
        printf("Failed to initialise window\n");
        free(app);
        return NULL;
    }

    app->renderer = SDL_CreateRenderer(app->window, -1, SDL_RENDERER_ACCELERATED);

    if(!app->renderer){
        printf("Failed to initialise renderer\n");
        free(app);
        return NULL;
    }

    app->texture = SDL_CreateTexture(app->renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, WINDOW_WIDTH, WINDOW_HEIGHT);

    if(!app->texture){
        printf("Failed to initialise buffer\n");
        free(app);
        return NULL;
    }

    app->buffer = malloc(sizeof(uint32_t) * WINDOW_HEIGHT * WINDOW_WIDTH);

    return app;
}

bool update_app(App *app, Fluid *fluid){

    // We calculate current FPS, from time elapsed, if not enough time elapsed since last time
    // image got refreshed we skip it to go calculate physics for more accurate results.
    float curr_fps = (float)SDL_GetPerformanceFrequency() / (SDL_GetPerformanceCounter() - app->last_tick); 
    if(curr_fps < app->fps){
        printf("FPS: %f\n", curr_fps);
        // get data from simulation to draw an image
        update_graphics(app, fluid);
        write_buffer_to_texture(app);
        SDL_RenderClear(app->renderer);
        SDL_RenderCopy(app->renderer, app->texture, NULL, NULL);
        SDL_RenderPresent(app->renderer);
        app->last_tick = SDL_GetPerformanceCounter();
    }


    // to be moved into handle_input function
    SDL_Event event;
    while (SDL_PollEvent(&event)){
        switch (event.type)
        {
             case SDL_QUIT:
                return false;
            default:
                break;
        }
    }

    return true;
}

void write_buffer_to_texture(App *app){
    uint32_t *pixel;
    int pitch;

    // In case that lock fails we omit copying data to texture (this could crash program idk).
    if(SDL_LockTexture(app->texture, NULL, (void **)&pixel, &pitch)) return;
    memcpy(pixel, app->buffer, sizeof(uint32_t) * WINDOW_HEIGHT * WINDOW_WIDTH);
    SDL_UnlockTexture(app->texture);
}

void update_graphics(App *app, Fluid *fluid){
    int chunk_width = WINDOW_WIDTH / fluid->width;
    int chunk_height = WINDOW_HEIGHT / fluid->height;

    for(int i = 0; i < fluid->width; i++){
        for(int j = 0; j < fluid->height; j++){
            for(int x = 0; x < chunk_width; x++){
                for(int y = 0; y < chunk_height; y++){
                    app->buffer[i * chunk_width + x + (y + j * chunk_height) * WINDOW_WIDTH] = 0xFF00 * fluid->density[i + j * fluid->width] / fluid->current_max_density + 0xFF;
                }
            }
        }
    }
    
}

void close_app(App *app){
    SDL_DestroyRenderer(app->renderer);
    SDL_DestroyWindow(app->window);
}
