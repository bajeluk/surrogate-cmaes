//#include "stdafx.h"
#include "mex.h" 
#include "math.h"
//#include "common.h"
//#include "RankSVMLearn.h"
//#include "dvector.h"


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

void normalizeX(double* x_Encoded, int ntrain, int nx, double* X_min, double* X_max)
{
	int i,j;
	double *xcur;
	for (j=0; j<nx; j++)
	{
		X_min[j] = x_Encoded[j];	
		X_max[j] = x_Encoded[j];
	}
	for (i=0; i<ntrain; i++)
	{
		xcur = &x_Encoded[i*nx];
		for(j=0; j<nx; j++)
		{
			if (xcur[j] < X_min[j])	X_min[j] = xcur[j];
			if (xcur[j] > X_max[j])	X_max[j] = xcur[j];
		}
	}
	for (i=0; i<ntrain; i++)
	{
		xcur = &x_Encoded[i*nx];
		for(j=0; j<nx; j++)
		{
			xcur[j]  = (xcur[j] - X_min[j] ) / (X_max[j] - X_min[j]);
		}
	}
}

void CalculateTrainingKernelMatrix(double* xvals, int npoints, int nx, double* K, 
								   double* TwoSigmaPow2, double sigma_A, double sigma_Pow, KernelType kernel, double kernelParam1, double kernelParam2)
{
	double distPow2,tmp;
	double *x1,*x2;
	double avrdist = 0;
	for (int i=0; i<npoints; i++)
		for (int j=i; j<npoints; j++)
		{
			if (kernel == kt_RBF)
			{
				if (i == j)	distPow2 = 0;
				else
				{
					x1 = &xvals[i*nx];
					x2 = &xvals[j*nx];
					distPow2 = 0;
					for (int k=0; k<nx; k++)
					{
						tmp = x1[k] - x2[k];
						distPow2 += tmp * tmp;
					}
				}
				K[i*npoints +j] = K[j*npoints +i] =distPow2;
				avrdist += sqrt(distPow2);
			}
			else
				K[i*npoints +j] = K[j*npoints +i] = MyKernel(kernel,&xvals[i*nx],&xvals[j*nx],nx,kernelParam1,kernelParam2);
		}

	avrdist /= ((npoints-1)*npoints/2);
	double sigma = sigma_A * pow(avrdist , sigma_Pow);
	*TwoSigmaPow2 = 2.0*sigma*sigma;

	double avrK = 0;
	for (int i=0; i<npoints; i++)
		for (int j=i; j<npoints; j++)
		{
			if (kernel == kt_RBF)
			{
				if (i == j)		K[i*npoints +j] = K[j*npoints +i] = 1;
				else			K[i*npoints +j] = K[j*npoints +i] = exp( -K[i*npoints +j] / *TwoSigmaPow2);
			}
		}
}



void OptimizeL(int ntrain, double* p_Ci, double epsilon, int niter,
						 double* p_Kij,double* p_dKij,double* p_alpha,double* p_sumAlphaDKij,double* p_div_dKij)
{
	int nAlpha = ntrain-1;
	double old_alpha, new_alpha, delta_alpha, sumAlpha, dL;
	int i,i1,j;

	double* p_dKii = new double[ntrain];

	for (i=0; i<nAlpha; i++)
	{
		for (j=0; j<nAlpha; j++)
		{
			p_dKij[i*nAlpha + j] = p_Kij[i*ntrain + j] - p_Kij[i*ntrain + (j+1)] - p_Kij[(i+1)*ntrain + j] + p_Kij[(i+1)*ntrain + (j+1)];
		}
		p_dKii[i] = p_dKij[i*nAlpha + i];
	}
	
	for (i=0; i<nAlpha;i++)
	{
	//	p_alpha[i] = p_Ci[i];// * rand()/(float)RAND_MAX;//p_Ci[i] * (0.95 + 0.05*rand()/(float)RAND_MAX);	// p_Ci[i] * rand()/(float)RAND_MAX;
	//	p_alpha[i] = p_Ci[i] * rand()/(float)RAND_MAX;
		p_alpha[i] = p_Ci[i] * (0.95 + 0.05*rand()/(float)RAND_MAX);
	}

	for (i=0; i<nAlpha; i++)
	{
		sumAlpha = 0;
		for (j=0; j<nAlpha;j++)
			sumAlpha += p_alpha[j] * p_dKij[i*nAlpha + j];
		p_sumAlphaDKij[i] = (epsilon - sumAlpha) / p_dKij[i*nAlpha + i];
	}

	for (i=0; i<nAlpha; i++)
		for (j=0; j<nAlpha; j++)
			p_div_dKij[i*nAlpha + j] = p_dKij[i*nAlpha + j] / p_dKij[j*nAlpha + j];


	double L = 0;
	double lastBestL = 0;
	double BestL = 0;

	int G = 0;
	int verb = 0;
	for (i=0; i<niter; i++)
	{	
		if (G == 0)
		{
			i1 = i%nAlpha;	//	int i1 = rand()%nAlpha;
			old_alpha = p_alpha[i1];
			new_alpha = old_alpha + p_sumAlphaDKij[i1];
			if (new_alpha > p_Ci[i1])		new_alpha = p_Ci[i1];
			if (new_alpha < 0)				new_alpha = 0;
			delta_alpha = new_alpha - old_alpha;


			//dL = delta_alpha * p_dKij[i1*nAlpha + i1] * ( p_sumAlphaDKij[i1] - 0.5*delta_alpha + 0);
			dL = delta_alpha * p_dKii[i1] * ( p_sumAlphaDKij[i1] - 0.5*delta_alpha + 0);

			if (dL > 0)
			{
				for (j=0; j<nAlpha; j++)
					p_sumAlphaDKij[j] -= delta_alpha * p_div_dKij[i1*nAlpha + j];

				p_alpha[i1] = new_alpha;
			}
			L = L + dL;
		}
		else
		{
			int jmax =0;
			double dLmax = -1e+30;
			for (int j1=0; j1<nAlpha; j1++)
			{
			
				i1 = j1;	//	int i1 = rand()%nAlpha;
				old_alpha = p_alpha[i1];
				new_alpha = old_alpha + p_sumAlphaDKij[i1];
				if (new_alpha > p_Ci[i1])		new_alpha = p_Ci[i1];
				if (new_alpha < 0)				new_alpha = 0;
				delta_alpha = new_alpha - old_alpha;


				//dL = delta_alpha * p_dKij[i1*nAlpha + i1] * ( p_sumAlphaDKij[i1] - 0.5*delta_alpha + 0);
				dL = delta_alpha * p_dKii[i1] * ( p_sumAlphaDKij[i1] - 0.5*delta_alpha + 0);
				if (dL > dLmax)
				{
					dLmax = dL;
					jmax = j1;
				}
			}
			i1 = jmax;	//	int i1 = rand()%nAlpha;
			old_alpha = p_alpha[i1];
			new_alpha = old_alpha + p_sumAlphaDKij[i1];
			if (new_alpha > p_Ci[i1])		new_alpha = p_Ci[i1];
			if (new_alpha < 0)				new_alpha = 0;
			delta_alpha = new_alpha - old_alpha;


			//dL = delta_alpha * p_dKij[i1*nAlpha + i1] * ( p_sumAlphaDKij[i1] - 0.5*delta_alpha + 0);
			dL = delta_alpha * p_dKii[i1] * ( p_sumAlphaDKij[i1] - 0.5*delta_alpha + 0);

			if (dL > 0)
			{
				for (j=0; j<nAlpha; j++)
					p_sumAlphaDKij[j] -= delta_alpha * p_div_dKij[i1*nAlpha + j];

				p_alpha[i1] = new_alpha;
			}
			L = L + dL;

			if (dL/(p_Ci[0]*p_Ci[0]) < 1e-20)
			{
			//	printf("%d %e\n",i,dL/(p_Ci[0]*p_Ci[0]));
								i = niter;
				verb = 0;
			}
		}


	}

	delete[] p_dKii;

}


void LearningRankN(double* x_training, double* x_trainingEncoded,   
		int nx, int ntrain, int niter, KernelType kernel, double epsilon, double* p_Ci, double* p_Cinv,
		double sigma_A, double sigma_Pow, double* xmean, int doEncoding, double* optAlphas, double* pTwoSigmaPow2, double kernelParam1, double kernelParam2,
		int normalize, double* X_min, double* X_max)

{
	
	//0. Init Temp Data
	int nAlpha = ntrain-1;	
	double* p_Kij = (double*)mxCalloc(ntrain*ntrain,sizeof(double));	double* p_dKij = (double*)mxCalloc(nAlpha*nAlpha,sizeof(double));
	double* p_alpha = (double*)mxCalloc(nAlpha,sizeof(double));		double* p_sumAlphaDKij = (double*)mxCalloc(nAlpha,sizeof(double));
	double* p_div_dKij = (double*)mxCalloc(nAlpha*nAlpha,sizeof(double));	double* p_dx = (double*)mxCalloc(nx,sizeof(double));
	double* Kvals = (double*)mxCalloc(ntrain,sizeof(double));
	int i,j;

	double ttotal = 0;
	//1. Transform training points to the new coordinate system and 
	//then compute Euclidean distance instead of 'EXPENSIVE' Mahalanobis distance
	if (doEncoding == 1)
	{	//p_Cinv is the C^-0.5 ;)
		Encoding(x_training, x_trainingEncoded, p_Cinv, xmean, ntrain, nx); 
	}

	if (normalize == 1)
	{
		normalizeX(x_trainingEncoded, ntrain, nx, X_min, X_max);
	}
	
	//2. Calculate the distance between points, then calculate sigma(gamma) and Kernel Matrix
	double TwoSigmaPow2 = 0;
	CalculateTrainingKernelMatrix(x_trainingEncoded, ntrain, nx, p_Kij, &TwoSigmaPow2, sigma_A, sigma_Pow, kernel, kernelParam1, kernelParam2);


	//3. Optimize alpha parameters
	OptimizeL(ntrain, p_Ci, epsilon, niter, p_Kij, p_dKij, p_alpha, p_sumAlphaDKij, p_div_dKij);

	
	*pTwoSigmaPow2 = TwoSigmaPow2;
	for (int i=0; i<nAlpha; i++)
	{
		optAlphas[i] = p_alpha[i];
	}

	mxFree(p_sumAlphaDKij);		mxFree(p_div_dKij);		mxFree(p_dKij);		mxFree(p_Kij);
	mxFree(p_dx);	mxFree(p_alpha);	mxFree(Kvals);
}





//[xtrainEncoded, alphas, TwoSigmaPow2] = RankSVMLearn(x_tr, nx, ntrain, niter, epsilon, Ci, kernel, ... 
//                                                        invsqrtC, sigmaA, sigmaPow, xmean, doEncoding, verbose);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//input
	int verbose = (int)(mxGetScalar(prhs[12]));

	double* X_training = mxGetPr(prhs[0]);
	int nx = (int)(mxGetScalar(prhs[1]));				//if (verbose == 1) printf("nx = %d\n",nx);
	int ntraining = (int)(mxGetScalar(prhs[2]));			//if (verbose == 1) printf("ntraining = %d\n",ntraining);
	int niter = (int)(mxGetScalar(prhs[3]));			//if (verbose == 1) printf("niter = %d\n",niter);
	double epsilon = (double)(mxGetScalar(prhs[4]));		//if (verbose == 1) printf("epsilon = %f\n",epsilon);
	double* p_Ci = mxGetPr(prhs[5]);
	int kernel = (int)(mxGetScalar(prhs[6]));			//if (verbose == 1) printf("kernel = %d\n",kernel);
	double* p_Cinv = mxGetPr(prhs[7]);
	double sigma_A = (double)(mxGetScalar(prhs[8]));		//if (verbose == 1) printf("sigma_A = %f\n",sigma_A);
	double sigma_Pow = (double)(mxGetScalar(prhs[9]));		//if (verbose == 1) printf("sigma_Pow = %f\n",sigma_Pow);
	double* xmean = mxGetPr(prhs[10]);
	int doEncoding = (int)(mxGetScalar(prhs[11]));			//if (verbose == 1) printf("doEncoding = %d\n",doEncoding);
	double KernelParam1 = (double)(mxGetScalar(prhs[12]));
	double KernelParam2 = (double)(mxGetScalar(prhs[13]));
	int normalize = (int)(mxGetScalar(prhs[14]));
	//
	int rowLen = mxGetN(prhs[0]);					//if (verbose == 1) printf("rowLen=%d\n",rowLen);
	int colLen = mxGetM(prhs[0]);					//if (verbose == 1) printf("colLen=%d\n",colLen);
	if ((ntraining != rowLen)||(nx != colLen)) 
	{	 // MatLab uses [i*colLen + j] notation, while we use [i*rowLen + j], so .. :) 
		//printf("Error: the matrix 'x_training' should have 'nx' rows and 'ntraining' columns");
		return;
	}

	//output
	plhs[0] = mxCreateDoubleMatrix(nx, ntraining, mxREAL);
	double* X_trainingEncoded = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(ntraining-1, 1, mxREAL);
	double* optAlphas = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* pTwoSigmaPow2 = mxGetPr(plhs[2]);

	plhs[3] = mxCreateDoubleMatrix(nx, 1, mxREAL);
	double* X_min = mxGetPr(plhs[3]);

	plhs[4] = mxCreateDoubleMatrix(nx, 1, mxREAL);
	double* X_max = mxGetPr(plhs[4]);

	LearningRankN(X_training, X_trainingEncoded, nx, ntraining, niter, (KernelType)kernel, epsilon, p_Ci, 
		p_Cinv, sigma_A, sigma_Pow, xmean, doEncoding, optAlphas, pTwoSigmaPow2, KernelParam1, KernelParam2, normalize, X_min, X_max);
}
