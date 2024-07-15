#include "vector.h"

vector2 vector_sum(vector2 v1, vector2 v2){
    vector2 ret;
    ret.x = v1.x + v2.x;
    ret.y = v1.y + v2.y;
    return ret;
}