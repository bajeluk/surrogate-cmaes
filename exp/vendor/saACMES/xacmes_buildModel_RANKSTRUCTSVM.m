function model = xacmes_buildModel_RANKSTRUCTSVM(cur_state,coeff)

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
            
epsilon = 1;
            
Cval = 10^6;    %default
Cval = 10^coeff(2);
nAlpha = nTrain - 1;
Ci = zeros(nAlpha,1);
z = 1:nAlpha;
           
rel = zeros(nTrain,nTrain);
Ci = zeros(nTrain,nTrain);
nConstr = 0;
nConstrToAdd = floor( nTrain*0 );
perc = nConstrToAdd / (nTrain*nTrain/2);
for i=1:nTrain
    for j=1:nTrain
%        if (j == i+1)   rel(i,j) = 1;   Ci(i,j) = Cval*(nTrain-i)^3.0; nConstr = nConstr + 1;
        if (j == i+1)   rel(i,j) = 1;   Ci(i,j) = Cval*(nTrain-i)^coeff(3); nConstr = nConstr + 1;
        else            rel(i,j) = 0;   Ci(i,j) = 0;  	end;
        if ((rand() < perc) && (j > i))
%            rel(i,j) = j-i;   Ci(i,j) = Cval*(j-i)^3.0; nConstr = nConstr + 1;
            rel(i,j) = j-i;   Ci(i,j) = Cval*(nTrain-i)^coeff(3); nConstr = nConstr + 1;
        end;
    end;
end;

niter = floor( 1000*nConstr );

if (nConstr > 5000)
    disp('LEARNING WILL FAIL IF WE CANNOT STORE nConstr X nConstr MATRIX IN MEMORY');
end;
%Ci(z) = Cval*((nAlpha-z).^3.0 );    %default
%Ci(z) = Cval*((nAlpha-z).^coeff(3) );
             
sigmaA = 1.0;   %default
sigmaA = coeff(4);
sigmaPow = 1.0;

doEncoding = 1.0;
verbose = 0;
x_tr = xtrain';
        
kernel = 0;
            
%[xtrainEncoded, alphas, TwoSigmaPow2] = RankSVMLearn(x_tr, N, nTrain, niter, epsilon, Ci, kernel, invsqrtC, sigmaA, sigmaPow, Xmean_model, doEncoding, verbose);
%
[xtrainEncoded, alphas, TwoSigmaPow2, nConstr] = RankStructSVMLearn(x_tr, rel, N, nTrain, niter, epsilon, Ci, kernel, ... 
                                                        invsqrtC, sigmaA, sigmaPow, Xmean_model, doEncoding, verbose);
 alphas = alphas(1:nConstr*3);

 model.modelType = 3;
 model.N = N;
 model.nTrain = nTrain;
 model.Xmean_model = Xmean_model;
 model.invsqrtC = invsqrtC;
 model.doEncoding = doEncoding;
 model.xtrainEncoded = xtrainEncoded;
 model.alphas = alphas;
 model.TwoSigmaPow2 = TwoSigmaPow2;
 model.nConstr = nConstr;
 
 model.nCrossValidation = nCrossValidation;
 model.CrossValidX = CrossValidX;
 model.CrossValidF = CrossValidF;
