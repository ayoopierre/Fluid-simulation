#ifndef MY_STENCIL_CPU
#define MY_STENCIL_CPU

#include <functional>

struct StencilCpu
{
    int offset;
    int di;
    int dj;
    std::function<double(int i, int j)> get_val;
};

typedef StencilCpu StencilCpu;

#endif