function model = xacmes_buildModel_SVR(cur_state,coeff)

global algo;

N = cur_state.N;
aSZ = algo.aSZ;
if (0)
    ARX = algo.ARX;
    ARF = algo.ARF;
end;
if (1)
    ARneval = algo.ARneval(1,1:aSZ);
    [ARneval, arindex] = sort(ARneval,2,'ascend');
    ARX = algo.ARX(:,arindex);
    ARF = algo.ARF(1,arindex);
end;
invsqrtC = cur_state.invsqrtC;
Xmean_model = cur_state.xmean;

nTrainMax = aSZ;
npts = floor(coeff(1));
if (nTrainMax > npts)
    nTrainMax = npts;
end;

[xtrain, ftrain, nTrain] = xacmes_selectTrainingPoints(ARX, ARF, aSZ, N, nTrainMax);    
if (nTrain < 0)
    disp( 'build model error: Ntrain < 0');
end;
            
nCrossValidation = 0;
CrossValidX = zeros(N,nCrossValidation);
CrossValidF = zeros(1,N);
for i=1:nCrossValidation
    index = 1 + floor( rand() * nTrain );
    nTrain = nTrain - 1;
    CrossValidX(:,i) = xtrain(index,:);
    CrossValidF(i) = ftrain(index);
    xtrain(index,:) = [];
    ftrain(index) = [];
end;
            
niter = floor( 10*nTrain);
%niter = floor( 10*nTrain*( 10^coeff(6)) );
epsilon = 1;

Cval = 10^6;    %default
Cval = 10^coeff(2);
nAlpha = nTrain;
Ci = zeros(nAlpha,1);
z = 1:nAlpha;
           
Ci(z) = Cval*((nAlpha-z).^3.0 );    %default
Ci(z) = Cval*((nAlpha-z).^coeff(3) );
             
sigmaA = 1.0;   %default
sigmaA = coeff(4);
sigmaPow = 1.0;

kernel = 0;
convEpsilon = 0;


doEncoding = 1.0;
verbose = 0;

x_tr = xtrain';

y_tr = ftrain;%
%y_tr = ftrain.^coeff(5);%
%y_tr = (1:nTrain).^coeff(5);%y_training;
y_tr = (1:nTrain);
%devEpsilon = (ftrain(2) - ftrain(1))*( 10^coeff(6));
devEpsilon = 0;

[xtrainEncoded, alphas, TwoSigmaPow2, bvalue] = SVRLearn(x_tr, y_tr, N, nTrain, niter, convEpsilon, devEpsilon, Ci, kernel, ... 
                                                        invsqrtC, sigmaA, sigmaPow, Xmean_model, doEncoding, verbose);

 model.modelType = 2;
 model.alphas = alphas;
 model.N = N;
 model.nTrain = nTrain;
 model.Xmean_model = Xmean_model;
 model.invsqrtC = invsqrtC;
 model.doEncoding = doEncoding;
 model.xtrainEncoded = xtrainEncoded;
 model.alphas = alphas;
 model.TwoSigmaPow2 = TwoSigmaPow2;
 model.bvalue = bvalue;
 
 model.nCrossValidation = nCrossValidation;
 model.CrossValidX = CrossValidX;
 model.CrossValidF = CrossValidF;
