#include "getResidual.h"

/*
  Calculate the residual with the U and F matrix
  INPUT: u, f, size
  OUTPUT: res
  size = N+1, size is the number of node in one direction

  More details see tutorial 10
 */


void getResidual(const double*u, const double*f, int size, double* res)
{
	int i, j;
	double h = 1.0/(size-1); // N = size-1
	for(i=1; i<size-1; i++)
		for(j=1; j<size-1; j++)
		{
			res[i*size+j] =  f[i*size+j] + (u[(i-1)*size+j] - 4*u[i*size+j] +  u[i*size + j-1] + u[i*size+j+1] + u[(i+1)*size+j])/(h*h);
		}
}
