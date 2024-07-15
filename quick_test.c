#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(){
    float *a = calloc(1, sizeof(float));
    int b = 1;
    *a = NAN;
    *a = *a;
    printf("%d", isnan(*a));
}