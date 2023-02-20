#include <stdio.h>
#include <math.h>

float distance(float A, float B, float C, float *a, float *b){

	float dx, dy, dz;
	float x,y,z;

	dx = fabsf(a[0] - b[0]);
	x = fmin(dx, fabsf(A - dx));
     
    dy = fabsf(a[1] - b[1]);
    y = fmin(dy, fabsf(B - dy));
     
    dz = fabsf(a[2] - b[2]);
    z = fmin(dz, fabsf(C - dz));

    return sqrt(x*x + y*y + z*z);
}