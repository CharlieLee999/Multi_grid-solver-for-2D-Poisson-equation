#include "gs.h"

void GS( const double *u0, const double *f, int size, int nIter, double *uNew)
{
	int iter, i, j;
	double h = 1.0/(size-1);
	double *uPrev = new double [size*size];
	double errMaxPos, errMinNeg, errMaxNorm;
	// Copy  all entries of u_0 to uNew as the inital value
	// Size = N+1 -> number of node in one direction
	for(i=0; i<size; i++)      
		for(j=0; j<size; j++)
		{
			uNew[i*size+j] = u0[i*size+j];
		}
	
	// GS-LEX
	for(iter=0; iter<nIter; iter++)   // Counter of iteration 
	{
		// Store the uNew in the last iteration as uPrev
		/*
		for(i=0; i<size*size; i++)		
		{
			uPrev[i] = uNew[i];
		}
		*/
		for(i=1; i< size-1; i++) // Iterate from i=1 to size-2(N-1) (total index-> from 0 to size-1)
		{
			for(j=1; j<size-1; j++)
			{
                // u[i][j] = 1/4( h*h* f[i][j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] )
                // Where u[i-1][j] and u[i][j-1] are already updated in previous iteration.
				uNew[i*size+j] =  1.0/4 * (h*h * f[i*size+j] + uNew[(i-1)*size+j] \
                       + uNew[i*size+j-1] + uNew[i*size+j+1] + uNew[(i+1)*size+j]);
			}
		}
		
		/*
		errMaxPos = 0;
		errMinNeg = 0;
		for(i=0; i<size*size; i++)
		{
			if(errMaxPos < uNew[i] - uPrev[i])
				errMaxPos = uNew[i] - uPrev[i];

			if(errMinNeg > uNew[i] - uPrev[i])
				errMinNeg = uNew[i] - uPrev[i];
		}

		if(abs(errMinNeg)> abs(errMaxPos))		
			errMaxNorm = abs(errMinNeg);	
		else 	
			errMaxNorm = abs(errMaxPos);
		*/
		
		// cout<< "Error L1 max norm of iteration "<< iter<<  " is "<< errMaxNorm<< endl;
		//cout<< "errMinNeg of iteration "<< iter<<  " is "<< errMinNeg<< endl;
		//cout<< "errMaxPos of iteration "<< iter<<  " is "<< errMaxPos<< endl;
	}
	//cout<< "Final error L1 max norm is "<< errMaxNorm<< endl;
	delete[] uPrev;
}
