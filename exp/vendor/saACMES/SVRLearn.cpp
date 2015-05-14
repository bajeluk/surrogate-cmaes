//#include "stdafx.h"

//#include "common.h"
#include "mex.h"
#include "math.h"
//#include "SVRLearn.h"


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


void LearningSVR(double* x_training, double* y_training_svm, double* x_trainingEncoded, 
		int nx, int ntrain, int niter, KernelType kernel, double convEpsilon, double devEpsilon, double* p_Ci, double* p_Cinv,
		double sigma_A, double sigma_Pow, double* xmean, int doEncoding, double* optAlphas, double* pTwoSigmaPow2, double* p_bvalue, double kernelParam1, double kernelParam2)

{

	double m_devEpsilon = devEpsilon;
//	double m_Cplus = 100;
//	double m_Cminus = 100;
	double m_epsilon = convEpsilon;

	int nImages = ntrain;
	double b_value = 0.0;
	double ret = 0.0;

	double* alpha = new double[2*nImages];
	double* linear = new double[2*nImages];
	double* lower = new double[2*nImages];
	double* upper = new double[2*nImages];
	double* gradient = new double[2*nImages];

	//	for (int i=0; i<ntrain; i++)
	//		y_training_svm[i] = -y_training_svm[i];
	//	Normalize(y_training_svm,ntrain);

	if (doEncoding == 1)
	{	//let's that we know that we want to do encoding so input invC is the C^-0.5 ;)
		Encoding(x_training, x_trainingEncoded, p_Cinv, xmean, ntrain, nx); 
	}

	double* K = new double[nImages * nImages];
	double TwoSigmaPow2 = 0;
	CalculateTrainingKernelMatrix(x_trainingEncoded, nImages, nx, K, &TwoSigmaPow2, sigma_A, sigma_Pow, kernel, kernelParam1, kernelParam2);

	int Dimension = 2*nImages;

	for(int i=0; i<nImages; i++)
	{
		alpha[i] = 0;
		alpha[i + nImages] = 0;;
		if (0)
		{
			alpha[i] = p_Ci[i]*rand()/(float)RAND_MAX;
			alpha[i + nImages] =  p_Ci[i]*rand()/(float)RAND_MAX;//m_Cplus-alpha[i];
		}		
		linear[i] = y_training_svm[i] - m_devEpsilon;
		linear[i + nImages] = y_training_svm[i] + m_devEpsilon;
		lower[i] = 0.0;
		lower[i + nImages] = -  p_Ci[i];
		upper[i] =  p_Ci[i];
		upper[i + nImages] = 0.0;
	}

	//Solve

	for (int i=0; i<2*nImages; i++)
		gradient[i] = linear[i];

	for (int i = 0; i < Dimension; i++)
	{
		double v = alpha[i];
		if (v != 0.0)
			for (int j = 0; j < Dimension; j++) 
				gradient[i] -= K[(i%nImages)*nImages+(j%nImages)] * v;
	}

	int iter = 0;
	int MaxIter = niter;

	while (iter != MaxIter)
	{
		int i1 = -1;
		int i2 = -1;
		int i1_prev,i2_prev;
		//MVP
		double largestUp = -1e100;
		double smallestDown = 1e100;
		for (int i = 0; i < Dimension; i++)
		{
			if (alpha[i] < upper[i])		
			{
				if (gradient[i] > largestUp)	// i такой, что Alpha(i) меньше Cmax(i) и имеем максимальный градиент
				{
					largestUp = gradient[i];
					if (i1 == -1)	i1_prev = nImages-1;
					else			i1_prev = i1;
					i1 = i;
				}
			}
			if (alpha[i] > lower[i])
			{
				if (gradient[i] < smallestDown)
				{
					smallestDown = gradient[i];	// j такой, что Alpha(j) больше Cmin(j) и имеем минимальный градиент
					if (i2 == -1)	i2_prev = nImages-1;
					else			i2_prev = i2;
					i2 = i;
				}
			}
		}

		if ( (i1%nImages) == (i2%nImages))
		{
			i1 = i1_prev;
		}
		int i1_mod = i1%nImages;
		int i2_mod = i2%nImages;

		double deltaGradient = largestUp - smallestDown;
		//	if (iter%1000 == 0)
		//		cout << iter << '\t' << deltaGradient << '\n';
		if (deltaGradient < m_epsilon)
			break;		// пока максимальный и минимальный градиент отлично более чем на epsilon

		// SMO update
		{
			double ai1 = alpha[i1];
			double ai2 = alpha[i2];
			double Li1 = lower[i1];
			double Ui1 = upper[i1];
			double Li2 = lower[i2];
			double Ui2 = upper[i2];

			// update alpha, that is, solve the sub-problem defined by i and j
			double nominator = gradient[i1] - gradient[i2];
			double denominator = K[i1_mod*nImages + i1_mod] + K[i2_mod*nImages + i2_mod] - 2.0 * K[i1_mod*nImages + i2_mod];
			//		double denominator = 2.0 - 2.0 * K[i1_mod*nImages + i2_mod];	//it works for RBF, because K(i,i) = 1 for RBF

			double mu = nominator / denominator;
			if (ai1 + mu < Li1) mu = Li1 - ai1;			//так если mu даст нам alpha < L, то ограничим его до L
			else if (ai1 + mu > Ui1) mu = Ui1 - ai1;	//так же ограничиваем сверху U
			if (ai2 - mu < Li2) mu = ai2 - Li2;			//тоже самое (по смыслу) для j
			else if (ai2 - mu > Ui2) mu = ai2 - Ui2;
			alpha[i1] += mu;
			alpha[i2] -= mu;

			// update the gradient
			double tmp;
			for (int i = 0; i < nImages; i++) 
			{
				tmp = mu * (K[i1_mod*nImages + i] - K[i2_mod*nImages + i]);
				gradient[i] -= tmp;
				gradient[nImages+i] -= tmp;
			}
		}

		iter++;
	}
	//	cout << "iIter:" << iter << '\n';

	//objective
	double objective = 0.0;
	for (int i = 0; i < Dimension; i++)
		objective += (gradient[i] + linear[i]) * alpha[i];
	objective =  0.5 * objective;

	// computation of b
	double value;
	double sum = 0.0;
	double lowerBound = -1e100;
	double upperBound = 1e100;
	unsigned int freeVars = 0;
	for (int i = 0; i < nImages; i++)
	{
		if (alpha[i] > 0.0)
		{
			value = gradient[i];
			if (alpha[i] <  p_Ci[i])
			{
				sum += value;
				freeVars++;
			}
			else
			{
				if (value > lowerBound) lowerBound = value;
			}
		}
		if (alpha[i + nImages] < 0.0)
		{
			value = gradient[i + nImages];
			if (alpha[i + nImages] > -  p_Ci[i])
			{
				sum += value;
				freeVars++;
			}
			else
			{
				if (value < upperBound) upperBound = value;
			}
		}
	}

	if (freeVars > 0)
		b_value = sum / freeVars;						// stabilized exact value
	else
		b_value = 0.5 * (lowerBound + upperBound);	// best estimate

	//alpha
	for (int i=0; i<nImages; i++)
		optAlphas[i] = alpha[i] + alpha[nImages+i];

	*p_bvalue = b_value;
	*pTwoSigmaPow2 = TwoSigmaPow2;

	delete[] K;
	delete[] lower;
	delete[] upper;
	delete[] gradient;
	delete[] linear;
	delete[] alpha;
}


//[xtrainEncoded, alphas, TwoSigmaPow2, bvalue] = SVRLearn(x_tr, y_tr, nx, ntrain, niter, convEpsilon, devEpsilon, Ci, kernel, ... 
//                                                        invsqrtC, sigmaA, sigmaPow, xmean, doEncoding, verbose);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//input
	int verbose = (int)(mxGetScalar(prhs[14]));

	double* X_training = mxGetPr(prhs[0]);
	double* Y_training = mxGetPr(prhs[1]);
	int nx = (int)(mxGetScalar(prhs[2]));				//if (verbose == 1) printf("nx = %d\n",nx);
	int ntraining = (int)(mxGetScalar(prhs[3]));			//if (verbose == 1) printf("ntraining = %d\n",ntraining);
	int niter = (int)(mxGetScalar(prhs[4]));			//if (verbose == 1) printf("niter = %d\n",niter);
	double convEpsilon = (double)(mxGetScalar(prhs[5]));		//if (verbose == 1) printf("convEpsilon = %f\n",convEpsilon);
	double devEpsilon = (double)(mxGetScalar(prhs[6]));		//if (verbose == 1) printf("devEpsilon = %f\n",devEpsilon);
	double* p_Ci = mxGetPr(prhs[7]);
	int kernel = (int)(mxGetScalar(prhs[8]));			//if (verbose == 1) printf("kernel = %d\n",kernel);
	double* p_Cinv = mxGetPr(prhs[9]);
	double sigma_A = (double)(mxGetScalar(prhs[10]));		//if (verbose == 1) printf("sigma_A = %f\n",sigma_A);
	double sigma_Pow = (double)(mxGetScalar(prhs[11]));		//if (verbose == 1) printf("sigma_Pow = %f\n",sigma_Pow);
	double* xmean = mxGetPr(prhs[12]);
	int doEncoding = (int)(mxGetScalar(prhs[13]));			//if (verbose == 1) printf("doEncoding = %d\n",doEncoding);
	double KernelParam1 = (double)(mxGetScalar(prhs[14]));
	double KernelParam2 = (double)(mxGetScalar(prhs[15]));
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

	plhs[1] = mxCreateDoubleMatrix(ntraining, 1, mxREAL);
	double* optAlphas = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* pTwoSigmaPow2 = mxGetPr(plhs[2]);

	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* pbvalue = mxGetPr(plhs[3]);

	LearningSVR(X_training, Y_training, X_trainingEncoded, nx, ntraining, niter, (KernelType)kernel, convEpsilon, devEpsilon, p_Ci, p_Cinv, sigma_A, sigma_Pow, xmean, doEncoding, optAlphas, pTwoSigmaPow2, pbvalue, KernelParam1, KernelParam2);
}
