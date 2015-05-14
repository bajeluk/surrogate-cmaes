function [xtrain, ftrain, nTrain] = xacmes_selectTrainingPoints(arrX, arrF, aSZ, xdim, nTrainMax)

nTrain = aSZ;
if (nTrain > nTrainMax)   nTrain = nTrainMax;   end;

indx1 = aSZ-nTrain+1:aSZ;
xtrain = arrX(:, indx1)';
ftrain = arrF(:, indx1)';

[ftrain, arindex] = sort(ftrain,'ascend');
xtrainnew = zeros(nTrain,xdim);
xtrainnew(:) = xtrain(arindex(1:nTrain),:);
xtrain = xtrainnew;

