#include "prolong.h"

/*
  This is the prolong function which prolongs a coarse mesh to a fine mesh by interpolating

  uC is the pointer of the coarse U matrix
  Nc is the number of piece of the coarse mesh, Nc = 2^{n-1} 
  uF is the pointer of the fine U matrix
  sizeC is the number of nodes of the coarse mesh, sizeC = Nc + 1
  sizeF is the number of nodes of the fine mesh,   sizeF = Nf + 1
 */


void PROLONG(const double* uC, int sizeC, double *uF)
{
	int i, j, ii, jj, idxI, idxJ, dist;
	int Nc = sizeC - 1;
	int Nf = 2 * Nc;
	int sizeF = Nf + 1;
	
	// Clear internal nodes of uF except the boundary node, because the boundary condition is zero. 
	for(i = 0; i < sizeF; i++)
		for(j = 0; j < sizeF; j++)
		{
			uF[i*sizeF + j] = 0;
		}
	
	// Update the uF with uC by scattering
	for(i=1; i<Nc; i++)		
	{
		ii = 2 * i;
		for(j=1; j<Nc; j++)
		{
			jj = 2 * j;
			// Calculate the coeffecient as -> pow(2, -dist)
			for(idxI=ii-1; idxI<ii+2; idxI++)
				for(idxJ=jj-1; idxJ<jj+2; idxJ++)
				{
					dist = abs(idxI-ii) + abs(idxJ-jj);
					uF[idxI*sizeF + idxJ] += pow(2, -dist) * uC[i*sizeC + j];
				}
		}
	}
}
