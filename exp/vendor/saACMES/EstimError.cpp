#include "mex.h" 
#include "math.h"


//err = EstimError(ranking, npoints);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//input
	double* ranking = mxGetPr(prhs[0]);
	int npoints = (int)(mxGetScalar(prhs[1]));

	//output
	plhs[0] = mxCreateDoubleMatrix(1,1,  mxREAL);

	
	double error = 0;
	for(int i=0; i<npoints; i++)
		for(int j=i+1; j<npoints; j++)
			if (ranking[i] >= ranking[j])
				error = error + 1;
	error = error / (npoints * (npoints - 1) / 2);
	double* err = mxGetPr(plhs[0]);
	*err = error;
}
