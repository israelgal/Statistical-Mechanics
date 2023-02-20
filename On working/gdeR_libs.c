#include <stdio.h>
#include <math.h>

//gcc -fPIC -shared -o gder_libs.so gder_libs.c

float distance(float A, float B, float C, float *a, float *b){
	
	float dx,dy,dz;
	float x,y,z;

	dx = fabsf(a[0] - b[0]);
    x = fmin(dx, fabsf(A - dx));
     
    dy = fabsf(a[1] - b[1]);
    y = fmin(dy, fabsf(B - dy));
     
    dz = fabsf(a[2] - b[2]);
    z = fmin(dz, fabsf(C - dz));

    return sqrt(x*x + y*y + z*z);
}




void RDF(float *rdf, float *coord, int amt_part, int A, int B, int C, int res, float dr){

	float dist;
	int index;



	for (int ii=0;ii<amt_part;ii++){
		for (int jj=ii;jj<amt_part;jj++){
			dist = distance(A,B,C,&coord[3*ii],&coord[3*jj]);
			index = (int)(dist/dr);
			if (0 < index && index < res) rdf[index] += 2.0;

		}
	}
}


