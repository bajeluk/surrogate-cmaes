
//#include "stdafx.h"
//#include "common.h"
#include "mex.h"
#include "math.h"
//#include "RankStructSVMLearn.h"
//#include "dvector.h"
//#include "RankStructSVMFunc.h"



struct Constr
{
	int i,j;
	double C;
	double epsilon;
	double CombXYXY;
	Constr()
	{}
};

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

void LearningRankStructSVM(double* x_training, double* x_trainingEncoded, double* rel,
				   int nx, int ntrain, int niter, KernelType kernel, double epsilon, double* p_Ci, double* p_Cinv,
				   double sigma_A, double sigma_Pow, double* xmean, int doEncoding, int* pnContrs, double* optAlphas, double* pTwoSigmaPow2, double kernelParam1, double kernelParam2)
{
	double* p_Kij = new double[ntrain*ntrain];

	if (doEncoding == 1)
	{	//p_Cinv = C^-0.5
		Encoding(x_training, x_trainingEncoded, p_Cinv, xmean, ntrain, nx); 
	}

	double TwoSigmaPow2 = 0.2;
	CalculateTrainingKernelMatrix(x_trainingEncoded, ntrain, nx, p_Kij, &TwoSigmaPow2, sigma_A, sigma_Pow,kernel, kernelParam1, kernelParam2);
	
	int nConstr = 0;
	for (int i=0; i<ntrain; i++)
		for (int j=0; j<ntrain; j++)
			if (rel[i*ntrain + j] != 0)	nConstr++;

	Constr* ActiveConstrs = new Constr[nConstr];
	double* p_alpha_constr = new double[nConstr];
	double* p_phi = new double[nConstr];
	double* pKconstr = new double[nConstr*nConstr];

	int iConstr = 0;
	for (int i=0; i<ntrain; i++)
		for (int j=0; j<ntrain; j++)
			if (rel[i*ntrain + j] != 0)
			{
				ActiveConstrs[iConstr].i = i;
				ActiveConstrs[iConstr].j = j;
				p_alpha_constr[iConstr] = p_Ci[i*ntrain + j];
				//p_alpha_constr[i] = C;//*rand()/(float)RAND_MAX;
				iConstr++;
			}

	//calculate phi
	for (int i=0; i<nConstr; i++)
	{
		int i1 = ActiveConstrs[i].i;
		int j1 = ActiveConstrs[i].j;
		double *pKpos = &pKconstr[i*nConstr];
		double phi = 0;
		for (int j=0; j<nConstr; j++)
		{
			pKpos[j] = CombKernel(p_Kij, ntrain, i1, j1, ActiveConstrs[j].i, ActiveConstrs[j].j);
		//	if (i == j)
		//	{
		//		double tt = p_Kij[i1*ntrain + j1] ;
		//		double delta = pKpos[j] - tt;
		//		double dt = delta - tt;
		//	}
			phi += p_alpha_constr[j] * pKpos[j];
		}
		p_phi[i] = phi;
	}


	double DL_sum = 0;
	double* pAlpha;

	for (int i=0; i<nConstr; i++)
	{
		Constr* pConstr = &ActiveConstrs[i];  
		int i1 = pConstr->i;
		int j1 = pConstr->j;
		pConstr->epsilon = rel[pConstr->i*ntrain + pConstr->j];
		pConstr->CombXYXY = pKconstr[i*nConstr + i];
		pConstr->C = p_Ci[i1*ntrain + j1];
	}

	for (int iter=0; iter<niter; iter++)
	{
		int iConstr = iter%nConstr;//rand()%nConstr;
		
		Constr* pConstr = &ActiveConstrs[iConstr];  
		double old_alpha = p_alpha_constr[iConstr];
		
		double epsilon_ij = pConstr->epsilon;//rel[pConstr->i*ntrain + pConstr->j];
		double CombXYXY =  pConstr->CombXYXY;//pKconstr[iConstr*nConstr + iConstr];
		double delta_alpha_opt = (epsilon_ij - p_phi[iConstr])/CombXYXY;
		double new_alpha = old_alpha + delta_alpha_opt;
		double C = pConstr->C;///p_Ci[i1*ntrain + j1];
		if (new_alpha > C)		new_alpha = C;
		if (new_alpha < 0)		new_alpha = 0;
		double delta_alpha = new_alpha - old_alpha;

		if (delta_alpha != 0)
		{
			double dL = delta_alpha * (epsilon_ij - p_phi[iConstr] - 0.5*delta_alpha*CombXYXY);
			if (dL > 0)
			{
				p_alpha_constr[iConstr] = new_alpha;
				pAlpha = &pKconstr[iConstr*nConstr];
				for (int i=0; i<nConstr; i++)
					p_phi[i] += delta_alpha*pAlpha[i];
			}
			DL_sum += dL;
		}
	}

	*pTwoSigmaPow2 = TwoSigmaPow2;
	*pnContrs = nConstr;
	for (int i=0; i<nConstr; i++)
	{
		optAlphas[3*i + 0] = ActiveConstrs[i].i;
		optAlphas[3*i + 1] = ActiveConstrs[i].j;
		optAlphas[3*i + 2] = p_alpha_constr[i];
	}

	delete []p_Kij;	delete []ActiveConstrs;	delete []p_alpha_constr;
	delete []p_phi;	delete []pKconstr;
}


double Lagrangian(double* alpha, double* p_Kij, double* rel, double* p_phi, int ntrain)
{
	double Lagr1 = 0;
	for (int i=0; i<ntrain; i++)
		for (int j=0; j<ntrain; j++)
			if (rel[i*ntrain + j] == 1)
				Lagr1 += alpha[i*ntrain + j];

	double Lagr2 = 0;
	for (int i=0; i<ntrain; i++)
		for (int j=0; j<ntrain; j++)
			if (rel[i*ntrain + j] == 1)
				for (int x=0; x<ntrain; x++)
					for (int y=0; y<ntrain; y++)
						if (rel[x*ntrain + y] == 1)
							Lagr2 += alpha[i*ntrain + j] * alpha[x*ntrain + y] * CombKernel(p_Kij,ntrain,i,j,x,y);
	double Lagr = Lagr1 - 0.5*Lagr2;
	return Lagr;
}

double Lagrangian2(double* p_alpha_constr,Constr* ActiveConstrs, int nConstr, double* pKconstr, int ntrain)
{
	double Lagr1 = 0;
	for (int i=0; i<nConstr; i++)
		Lagr1 += p_alpha_constr[i];

	double Lagr2 = 0;
	for (int i=0; i<nConstr; i++)
		for (int j=0; j<nConstr; j++)
			Lagr2 += p_alpha_constr[i] * p_alpha_constr[j] * pKconstr[i*nConstr + j];

	double Lagr = Lagr1 - 0.5*Lagr2;
	return Lagr;
}

//[xtrainEncoded, alphas, TwoSigmaPow2, nConstr] = RankStructSVMLearn(x_tr, rel, nx, ntrain, niter, epsilon, Ci, kernel, ... 
//                                                        invsqrtC, sigmaA, sigmaPow, xmean, doEncoding, verbose);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//input
	int verbose = (int)(mxGetScalar(prhs[13]));

	double* X_training = mxGetPr(prhs[0]);
	double* rel = mxGetPr(prhs[1]);
	int nx = (int)(mxGetScalar(prhs[2]));				//if (verbose == 1) printf("nx = %d\n",nx);
	int ntraining = (int)(mxGetScalar(prhs[3]));			//if (verbose == 1) printf("ntraining = %d\n",ntraining);
	int niter = (int)(mxGetScalar(prhs[4]));			//if (verbose == 1) printf("niter = %d\n",niter);
	double epsilon = (double)(mxGetScalar(prhs[5]));		//if (verbose == 1) printf("epsilon = %f\n",epsilon);
	double* p_Ci = mxGetPr(prhs[6]);
	int kernel = (int)(mxGetScalar(prhs[7]));			//if (verbose == 1) printf("kernel = %d\n",kernel);
	double* p_Cinv = mxGetPr(prhs[8]);
	double sigma_A = (double)(mxGetScalar(prhs[9]));		//if (verbose == 1) printf("sigma_A = %f\n",sigma_A);
	double sigma_Pow = (double)(mxGetScalar(prhs[10]));		//if (verbose == 1) printf("sigma_Pow = %f\n",sigma_Pow);
	double* xmean = mxGetPr(prhs[11]);
	int doEncoding = (int)(mxGetScalar(prhs[12]));			//if (verbose == 1) printf("doEncoding = %d\n",doEncoding);
	double KernelParam1 = (double)(mxGetScalar(prhs[13]));
	double KernelParam2 = (double)(mxGetScalar(prhs[14]));

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

	plhs[1] = mxCreateDoubleMatrix(ntraining*ntraining*3, 1, mxREAL);
	double* optAlphas = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* pTwoSigmaPow2 = mxGetPr(plhs[2]);

	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* pnConstr = mxGetPr(plhs[3]);
	int nConstr = 0;

	LearningRankStructSVM(X_training, X_trainingEncoded, rel, nx, ntraining, niter, (KernelType)kernel, epsilon, p_Ci, p_Cinv, sigma_A, sigma_Pow, xmean, doEncoding, &nConstr, optAlphas, pTwoSigmaPow2, KernelParam1, KernelParam2);
	*pnConstr = nConstr;
}

