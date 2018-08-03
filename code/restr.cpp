#include "restr.h"
/*
	uF(input) -> fine mesh, [Nf + 1] * [Nf + 1]
	uC(output) -> Coarse mesh, [Nc + 1] * [Nc + 1]
	In this function, we only calculate the interal nodes of coarse mesh because the nodes on
Boundary will keep fixed as zero.
	Use the algorithm in Tutorial 10
*/
void RESTR(const double* uF, int sizeF, double *uC)
{	
	int i, j, ii, jj;
	int Nf = sizeF - 1;
	int Nc = Nf/2;
	int sizeC = Nc + 1;
	// Nc = 2^{n-1}
	// Nf = 2^{n}
	
	/*	Only update the interal nodes of the coarse mesh.
		Interal nodes -> [1, 2, ..., Nc-1] number = Nc-1
		Total nodes   -> [0, 1, ..., Nc] number = Nc+1
	*/
	for(i=1; i<Nc; i++)
	{
		ii = 2 * i;
		for(j=1; j<Nc; j++)
		{
			jj = 2 * j;
            
			uC[i*sizeC + j] = 1.0/16 * ( uF[(ii-1)*sizeF + jj-1] + 2*uF[ii*sizeF + jj-1]
			+ uF[(ii+1)*sizeF + jj-1] + 2*uF[(ii-1)*sizeF + jj] + 4*uF[ii*sizeF + jj]
			+ 2*uF[(ii+1)*sizeF + jj] + uF[(ii-1)*sizeF + jj+1] + 2*uF[ii*sizeF + jj+1]
			+ uF[(ii+1)*sizeF + jj+1] );
		}
	}
}
