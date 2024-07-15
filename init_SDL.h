#include <stdbool.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <string.h>
#include "fluid.h"

#define WINDOW_HEIGHT 600
#define WINDOW_WIDTH 800
#define FPS 60

//TODO -> change name to smth like handle_SDL.h

typedef struct{
    SDL_Renderer *renderer;
    SDL_Window *window;
    SDL_Texture *texture;
    uint32_t *buffer;
    Uint64 last_tick;
    int last_x_mouse;
    int last_y_mouse;
    int fps;
} App;

App *init_app();
bool update_app(App *app, Fluid *fluid);
void close_app(App *app);
void write_buffer_to_texture(App *app);
void update_graphics(App *app, Fluid *fluid);
