#include "globalHeader.h"
#include "gs.h"
#include "restr.h"
#include "getResidual.h"
#include "prolong.h"
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
void mg(const double* uIn, const double *f, int l, int gama, int nu1, int nu2, int sizeF, double *uOut);
void checkGS();
void checkRESTR_PROLONG();
// Please compile with the flag -std=c++11
// g++ -std=c++11 *.cpp -o test


int main( int argc, char*argv[])
{
    // checkGS();
    // checkRESTR_PROLONG();
    
    
    // parameters 
    int n = atoi(argv[1]);      // n -> N = 2^{n} 
    int nIter = 20; // total number of iterations
    int l = n-1;
    int gamma = 2;
    int nu1 = atoi(argv[2]); // 1 or 2 or any positive integer.
    int nu2 = atoi(argv[3]);
    /* l must be n-1 to use as many multigrid as possible
       and then the coarest mesh has only one interior node.
     */

    // csv File writer
    ofstream outfile;
	int Nf, Nc, sizeF, sizeC;
    Nf = pow(2, n);
    Nc = Nf/2;
    sizeF = Nf + 1;
    sizeC = Nc + 1;
    double h = 1.0/Nf;
    int i, j, iter;
    double resInit, errInit;
    double *resNormVec, *resRatio;
    double *uIn, *uFExact, *f, *uOut, *resIter;
    double *errNormVec, *errRatio, *errIter;
    
    
    uIn = new double[sizeF*sizeF](); // Initialization
    uFExact = new double[sizeF*sizeF];
    uOut = new double[sizeF*sizeF](); // Initialization
    f = new double[sizeF*sizeF];
    resNormVec = new double[nIter+1];
    resRatio = new double[nIter+1];
    resIter = new double[sizeF*sizeF];

    errNormVec = new double[nIter+1];
    errRatio = new double[nIter+1];
    errIter = new double[sizeF*sizeF];
    
    // Calculate the exact fine mesh
    for(i=0; i<sizeF; i++)
        for(j=0; j<sizeF; j++)
        {
            /*
              uExact = sin(2*pi*x) * sin( 2*pi*y)
            */
            uFExact[i*sizeF+j] = sin(2*M_PI* h*i) * sin(2*M_PI* h*j);
            f[i*sizeF+j] = 8*M_PI*M_PI * uFExact[i*sizeF+j];
        }

    // initial residual res0 = f because u0 = 0;
    resInit = *max_element(f, f+sizeF*sizeF);
    resNormVec[0] = resInit;
    resRatio[0] = 1;

    // initial error err0 = uFExact because u0 = 0    
    errInit = *max_element(uFExact, uFExact+sizeF*sizeF);
    errNormVec[0] = errInit;
    errRatio[0] = 1;
    
    // Multi-grid iteration
    for(iter=0; iter<nIter; iter++)
    {
    	cout<< "\t\t Iteration-> -> " << iter<< endl;
        mg(uIn, f, l, gamma, nu1, nu2, sizeF, uOut);

        // Calculate error matrix and  the max norm of error
        for(i=0; i<sizeF*sizeF; i++)
        {
            errIter[i] = abs(uOut[i] - uFExact[i]);
        }
        errNormVec[iter+1] = *max_element(errIter, errIter+sizeF*sizeF);
        errRatio[iter+1] = errNormVec[iter+1]/errInit;
        
        // Calculate the residual matrix and the max norm of residual;
        getResidual(uOut, f, sizeF, resIter);
        for(i=0; i< sizeF*sizeF; i++)
        {
            resIter[i] = abs(resIter[i]);
        }
        resNormVec[iter+1] = *max_element(resIter, resIter+sizeF*sizeF);
        resRatio[iter+1] = resNormVec[iter+1] / resInit;
        cout<< "resRatio[" << iter+1<< "] is "<< resRatio[iter+1]<< endl;
        
        // Copy the uOut to uIn and begin preceding iteration.
        for(i=0; i<sizeF*sizeF; i++)
        {
            uIn[i] = uOut[i];
        }
    }

    
    cout<< "Writing to csv" <<endl;
    string fNameErr = "result/errNormVec_n_" + to_string(n)+"_nu1_"+to_string(nu1)+"_nu2_"+to_string(nu2)+".csv";
    outfile.open(fNameErr);
    for(i=0; i<nIter+1; i++)
    {
        outfile<< errNormVec[i]<< ","<< i<< "," <<"\n";
    }
    outfile.close();

    string fNameErrRatio = "result/errRatio_n_" + to_string(n)+"_nu1_"+to_string(nu1)+"_nu2_"+to_string(nu2)+".csv";
    outfile.open(fNameErrRatio);
    for(i=0; i<nIter+1; i++)
    {
        outfile<< errRatio[i]<< ","<< i<< ","  <<"\n";
    }
    outfile.close();
    
    string fNameRes = "result/resNormVec_n_" + to_string(n)+"_nu1_"+to_string(nu1)+"_nu2_"+to_string(nu2)+ ".csv";
    outfile.open(fNameRes);
    for(i=0; i<nIter+1; i++)
    {
        outfile<< resNormVec[i]<< ","<< i<< ","  <<"\n";
    }
    outfile.close();

    string fNameResRatio = "result/resRatio_n_" + to_string(n)+"_nu1_"+to_string(nu1)+"_nu2_"+to_string(nu2)+".csv";
    outfile.open(fNameResRatio);
    for(i=0; i<nIter+1; i++)
    {
        outfile<< resRatio[i]<< "," << i<< "," <<"\n";
    }
    outfile.close();
    
    delete[] resNormVec;
    delete[] resRatio;
    delete[] uIn;
    delete[] uFExact;
    delete[] uOut;
    delete[] f;
    delete[] resIter;
    
	return 0;
}


void mg(const double* uIn, const double *f, int l, int gama, int nu1, int nu2, int sizeF, double *uOut)
{
	//cout<< "Calling MG function"<< endl;
	double *uSmoothed;
	double *resF, *resC, *uPlus, *errTildeF;
	double *errTildeC, *errTmp;

	int i, j;

	int sizeC = (sizeF-1)/2 + 1; // Nf = sizeF - 1, Nc = Nf/2, sizeC = Nc + 1

	uSmoothed = new double[sizeF*sizeF];
	resF = new double[sizeF*sizeF];
	resC = new double[sizeC*sizeC];
	uPlus = new double[sizeF*sizeF];
	errTildeF = new double[sizeF*sizeF];
	errTildeC = new double[sizeC*sizeC]();
	errTmp = new double[sizeC*sizeC]();

	/*
	if(l==6) // n-1 The first step of the MG function
	{
		cout<< "Current resC pointer is "<< resC<< endl;
	}
	*/
	//cout<< "Calling GS function"<< endl;
	GS(uIn, f, sizeF, nu1, uSmoothed);

	//cout<< "Calling getResidual function" << endl;
	getResidual(uSmoothed, f, sizeF, resF);

	//cout<< "Calling the restriction function" << endl;
	RESTR(resF, sizeF, resC);

	// Let coarse residual change the sign because -resC is the argument of GS and GM functions.
	for(i=0; i<sizeC*sizeC; i++)
	{
		resC[i] = -resC[i];
	}

	if(l==1)
	{
	/*
	Condition: l=n-1, sizeF = 2^{n} + 1.
	If l=n-1, the coarest mesh only has one interior node.
	Why l = n-1:
	Nf = 2^{n},
	and each restriction reduce N by half, then Nc = 2^{n-1}.
	After l-1 iteration, the coarest Nc is 2^{1},
	then there are only two boundary nodes and one interior node.
	Now, we can solve the error function exactly with the one-step
	Gauss-seidel funtion.
	*/
	//  cout<< "\n\tExactly solving the error, sizeC ->  "<< sizeC<< endl;
	//errTildeC[0] = 0; // scalar
		GS(errTildeC, resC, sizeC, 1, errTmp);
		for(i=0; i<sizeC*sizeC; i++)
		{
			errTildeC[i] = errTmp[i];
		}
	}
	else
	{
		for(j=0; j<gama; j++)
		{
			mg(errTildeC, resC, l-1, gama, nu1, nu2, sizeC, errTmp);

			// Copy the errTmp from the mg function to errTildeC
			for(i=0; i< sizeC*sizeC; i++)
			{
				errTildeC[i] = errTmp[i];
			}
		}
	}

	//cout<< "Calling prolongation function,  l-> " << l<<endl;
	PROLONG(errTildeC, sizeC, errTildeF);

	//cout<< "Calculating uPlus"<< endl;
	for(i=0; i<sizeF*sizeF; i++)
	{
		uPlus[i] = uSmoothed[i] - errTildeF[i];        
	}

	//cout<< "Calling GS function" << endl;
	GS(uPlus, f, sizeF, nu2, uOut);

	/*
	if(l==3) // n-1 The final step of MG function
	{
		cout<< "Current resC pointer is "<< resC<< endl;        
	}
	*/
	//cout<< "Clearing memory,  l-> " << l<<endl;
	delete[] uSmoothed;
	// cout<< "Clearing memory resC,  l-> " << l<<endl;
	delete[] resC;
	// cout<< "Clearing memory errTildeC,  l-> " << l<<endl;
	delete[] errTildeC;
	// cout<< "Clearing memory errTildeF,  l-> " << l<<endl;
	delete[] errTildeF;
	// cout<< "Clearing memory errTmp,  l-> " << l<<endl;
	delete[] errTmp;
	// cout<< "Clearing memory uPlus,  l-> " << l<<endl;
	delete[] uPlus;
	// cout<< "Clearing memory resF,  l-> " << l<<endl;
	delete[] resF;
}



void  checkGS()
{
    double *uInit, *uExact, *f, *uSmoothed, *uDiff;
    int N = 10;
    int sizeF = N + 1;
    double h = 1.0/N;
    int i, j;

    // Calculate the error 
    double errMaxPos, errMinNeg, errMaxNorm;
    
    uInit = new double[sizeF*sizeF]();
    uExact = new double[sizeF*sizeF];
    f = new double[sizeF*sizeF];
    uSmoothed = new double[sizeF*sizeF]();
    uDiff = new double[sizeF*sizeF];
    
    for(i=0; i<sizeF; i++)
        for(j=0; j<sizeF; j++)
        {
            /*
              uExact = sin(2*pi*x) * sin( 2*pi*y)
              f = 8*pi*pi * sin(2*pi*x) * sin(2*pi*y)
                = 8*pi*pi * uExact
             */
            uExact[i*sizeF+j] = sin(2*M_PI* h*i) * sin(2*M_PI* h*j);
            f[i*sizeF+j] = 8*M_PI*M_PI * uExact[i*sizeF+j];
        }
    cout<< setprecision(9)<< "M_PI is "<< M_PI<< endl;
    
    cout<< "uExact"<< endl;
    /*
    for(i=0; i<sizeF; i++)
	{
		for(j=0; j<sizeF; j++)
        {
        	cout<< "["<< i<< "]["<< j<< "] is "<<uExact[i*sizeF+j] << "  ";
        }
     	cout<< "\n" << endl;
    }
    */
    // N = 10 -> nIter = 151
    // N = 100 -> nIter = 5000+
    GS(uInit, f, sizeF, 151, uSmoothed);
    
    for(i=0; i<sizeF*sizeF; i++)
    {
    	uDiff[i] = uSmoothed[i] - uExact[i];
    }
    
    errMaxPos = *max_element(uDiff, uDiff + sizeF*sizeF);
    errMinNeg = *min_element(uDiff, uDiff + sizeF*sizeF);
    
    if(abs(errMinNeg)> abs(errMaxPos))		
			errMaxNorm = abs(errMinNeg);	
		else 	
			errMaxNorm = abs(errMaxPos);
	
	cout<< "Converged errMaxPos between uSmoothed and uExact is" << errMaxPos<< endl;
	cout<< "Converged errMinNeg between uSmoothed and uExact is" << errMinNeg<< endl;
    cout<< "Converged max error between uSmoothed and uExact is" << errMaxNorm<< endl;
    
    delete[] uInit;
    delete[] uExact;
    delete[] f;
    delete[] uSmoothed;
    delete[] uDiff;
}

void checkRESTR_PROLONG()
{
	int n = 7;
	int Nf, Nc, sizeF, sizeC;
    Nf = pow(2, n);
    Nc = Nf/2;
    sizeF = Nf + 1;
    sizeC = Nc + 1;
    double h = 1.0/Nf;
    double H = 1.0/Nc;
    int i, j;

    double *uF, *uC, *uCExact, *uCDiff, *uFExact, *uFDiff;

    uF = new double[sizeF*sizeF](); // Initialization
    uFExact = new double[sizeF*sizeF];
    uFDiff = new double[sizeF*sizeF];
    
    uC = new double[sizeC*sizeC](); // Initialization
    uCExact = new double[sizeC*sizeC];
    uCDiff = new double[sizeC*sizeC];

    // Calculate the exact fine mesh
    for(i=0; i<sizeF; i++)
        for(j=0; j<sizeF; j++)
        {
            /*
              uExact = sin(2*pi*x) * sin( 2*pi*y)
             */
            uFExact[i*sizeF+j] = sin(2*M_PI* h*i) * sin(2*M_PI* h*j);
        }
    // Calculate the exact coarse mesh
    for(i=0; i<sizeC; i++)
        for(j=0; j<sizeC; j++)
        {
            uCExact[i*sizeC+j] = sin(2*M_PI *H*i) * sin(2*M_PI * H*j);
        }

    // RESTRICTION OPREATOR
    RESTR(uFExact, sizeF, uC);
    
    for(i=0; i<sizeC*sizeC; i++)
    {
        uCDiff[i] = abs(uC[i] - uCExact[i]);
    }
    
    cout<< "Max error between the restricted and exact coarse mesh is "<< *max_element(uCDiff, uCDiff+sizeC*sizeC)<< endl;

    // PROLONG OPERATOR
    PROLONG(uCExact, sizeC, uF);
    for(i=0; i<sizeF*sizeF; i++)
    {
        uFDiff[i] = abs(uF[i] - uFExact[i]);
    }
    cout<< "Max error between the prolonged and exact fine mesh is "<< *max_element(uFDiff, uFDiff+sizeF*sizeF)<< endl;
    
    delete[] uF;
    delete[] uFExact;
    delete[] uFDiff;
    delete[] uC;
    delete[] uCDiff;
    delete[] uCExact;
}
