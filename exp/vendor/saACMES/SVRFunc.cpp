//#include "stdafx.h"
#include "mex.h" 
#include "math.h"
//#include "common.h"
//#include "SVRFunc.h"

enum KernelType
{
	kt_RBF = 0,
	kt_Polynomial = 1,
};

double DistPow2_Euclidean(double* x1, double* x2,int nx)
{
	double tmp;
	double dist = 0;
	for (int i=0; i<nx; i++)
	{
		tmp = x1[i] - x2[i];
		dist += tmp * tmp;
	}
	return dist;
}

double Kernel_RBF(double* x1, double* x2,int nx, double TwoSigmaPow2)
{
	return exp(- DistPow2_Euclidean(x1, x2, nx) / TwoSigmaPow2);
}

double Kernel_Polynomial(double* x1, double* x2,int nx, double c, double d)
{
	double KK = 0;
	for (int k=0; k<nx; k++)
		KK += x1[k] * x2[k];
	return pow(KK + c,d);
}

double MyKernel(KernelType kernel, double* x1, double* x2,int nx, double param1, double param2)
{
	if (kernel == kt_RBF)			return Kernel_RBF(x1,x2,nx,param1);
	if (kernel == kt_Polynomial)		return Kernel_Polynomial(x1,x2,nx,param1,param2);
	return -1;
}

double CombKernel(double* p_Kij, int ntrain, int i, int j, int x, int y)
{
	return (p_Kij[i*ntrain + x] - p_Kij[i*ntrain + y] - p_Kij[j*ntrain + x] + p_Kij[j*ntrain + y]);
}


void Encoding(double* x, double* x_Encoded, double* invsqrtC, double* xmean, int npoints, int nx)	//x'(i) = C^(-0.5) * ( x(i) - xmean(i) )
{
	double *xcur, *xnew;
	double sum;
	int jrow;
	double* dx = (double*)mxCalloc(nx,sizeof(double));
	for (int i=0; i<npoints; i++)
	{
		xcur = &x[i*nx];
		xnew = &x_Encoded[i*nx];
		for (int j=0; j<nx; j++)
			dx[j] = xcur[j] - xmean[j];
		for (int j=0; j<nx; j++)
		{
			sum = 0;
			jrow = j*nx;
			for (int k=0; k<nx; k++)
				sum += invsqrtC[jrow + k] * dx[k];
			xnew[j] = sum;
		}
	}
	mxFree(dx);
}

void EvaluatePointSVR(double* Fit, double* x_training, double* x_test, int ntest, int nx, int ntrain, double* p_alpha, double TwoSigmaPow2, double* p_Cinv, double* xmean, int doEncoding, double b_value, KernelType kernel, double kernelParam1, double kernelParam2)
{
	double* x_Encoded = (double*)mxCalloc(nx,sizeof(double));

	for(int i=0; i<ntest; i++)
	{
		double* curX = &x_test[i*nx];
		if (doEncoding == 1)
		{	
			Encoding(&x_test[i*nx], x_Encoded, p_Cinv, xmean, 1, nx); 
			curX = x_Encoded;
		}

		double distPow2;	
		double curFit = b_value;
		for (int j=0; j<ntrain; j++)
			if (p_alpha[j] != 0)
			{
				//distPow2 = DistPow2_Euclidean(curX, &x_training[j*nx], nx);
				//curFit += p_alpha[j] * exp(- distPow2 / TwoSigmaPow2);
				curFit += p_alpha[j] * MyKernel(kernel,curX, &x_training[j*nx], nx, kernelParam1,kernelParam2);
			}
		Fit[i] = curFit;
	}	
	mxFree(x_Encoded);
}


//y = SVRFunc(model.xtrainEncoded, x_test, npoints, model.N, ... 
//          	model.nTrain, model.alphas, model.TwoSigmaPow2, model.invsqrtC, model.Xmean_model, model.doEncoding, model.bvalue);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//input
	double* X_training = mxGetPr(prhs[0]);
	double* X_test = mxGetPr(prhs[1]);
	int ntest = (int)(mxGetScalar(prhs[2]));
	int nx = (int)(mxGetScalar(prhs[3]));				
	int ntraining = (int)(mxGetScalar(prhs[4]));			
	double* p_alpha = mxGetPr(prhs[5]);
	double TwoSigmaPow2 = (double)(mxGetScalar(prhs[6]));	
	double* p_Cinv = mxGetPr(prhs[7]);
	double* xmean = mxGetPr(prhs[8]);
	int doEncoding = (int)(mxGetScalar(prhs[9]));
	double b_value = (double)(mxGetScalar(prhs[10]));
	int kernel = (int)(mxGetScalar(prhs[11]));
	KernelType Kernel = (KernelType)kernel;
	double KernelParam1 = (double)(mxGetScalar(prhs[12]));
	double KernelParam2 = (double)(mxGetScalar(prhs[13]));


	int rowLen = mxGetN(prhs[0]);						
	int colLen = mxGetM(prhs[0]);							
	if ((ntraining != rowLen)||(nx != colLen)) 
	{	 // MatLab uses [i*colLen + j] notation, while we use [i*rowLen + j], so .. :) 
	//printf("Error: the matrix 'x_training' should have 'nx' rows and 'ntraining' columns");
	return;
	}
	//output
	plhs[0] = mxCreateDoubleMatrix(1,ntest,  mxREAL);
	double* Fit = mxGetPr(plhs[0]);

	EvaluatePointSVR(Fit, X_training, X_test, ntest,nx, ntraining, p_alpha, TwoSigmaPow2, p_Cinv, xmean, doEncoding, b_value,Kernel, KernelParam1, KernelParam2);
}
