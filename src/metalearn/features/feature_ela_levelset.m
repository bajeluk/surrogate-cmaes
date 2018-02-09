function ft = feature_ela_levelset(X, y)
% ft = FEATURE_ELA_LEVELSET(X, y) returns ELA levelset features
% for dataset [X, y]. 
% 
% The dataset [X, y] is split into two classes by a specific treshold
% calculated using different quantiles. Linear and quadratic discriminant 
% analysis are used to predict whether the objective values y fall below 
% or exceed the calculated threshold. The extracted features are based on 
% the distribution of the resulting cross-validated mean misclassification 
% errors of each classifier.
%
% Quantiles: 0.1, 0.25, 0.5
% Methods:   lda, qda (future work: mda)
% 
% Features:
%   mmce_[method]_[quantile]       - mean misclassification error of 
%                                    appropriate method using defined 
%                                    quantile
%   [method1]_[method2]_[quantile] - ratio between mean misclassification 
%                                    errors of two methods using defined 
%                                    quantile

  if nargin < 2
    help feature_ela_levelset
    if nargout > 0
      ft = struct();
    end
    return
  end
  
  % initialize
  qnt = [0.1, 0.25, 0.5];
  classMethods = {'lda', 'qda'};
  nFolds = 10;
  
  nData = numel(y);
  nQuant = numel(qnt);
  nClassMet = numel(classMethods);
  
  % check number of folds
  if nData < nFolds
    nFolds = nData;
  end
  % create instances for 10-fold CV
  cp = cvpartition(nData, 'KFold', nFolds);
  % quantile tresholds
  y_tresh = quantile(y, qnt);
  mmce = zeros(nQuant, nClassMet);
  % quantile loop
  for q = 1:nQuant
    % class labels according to quantiles
    y_class = (y < y_tresh(q));
    for cm = 1:nClassMet
      switch classMethods{cm}
        case 'lda'
          classFcn = @(xTrain, yTrain, xTest) oneFold(xTrain, yTrain, xTest, 'linear');
        case 'qda'
          classFcn = @(xTrain, yTrain, xTest) oneFold(xTrain, yTrain, xTest, 'quadratic');
        % TODO: mda
        otherwise
          error('No such classification method as %s available for ELA levelset features implemented', ...
                classMethods{cm})
      end
      % CV
      mmce(q, cm) = crossval('mcr', X, y_class, 'predfun', classFcn, 'partition', cp);
      % mean misclassification error  features
      mmceName = sprintf('mmce_%s_%d', classMethods{cm}, 100*qnt(q));
      ft.(mmceName) = mmce(q, cm);
    end
    % mmce combination features
    combId = nchoosek(1:nClassMet, 2);
    for ci = 1:size(combId, 1)
      combiName = sprintf('%s_%s_%d', ...
                          classMethods{combId(ci, 1)}, ...
                          classMethods{combId(ci, 2)}, ...
                          100*qnt(q));
      ft.(combiName) = mmce(q, combId(ci, 1)) / mmce(q, combId(ci, 2));
    end
  end
  
end

function yTest = oneFold(xTrain, yTrain, xTest, type)
% one fold function
  try
    mdl = fitcdiscr(xTrain, yTrain, 'DiscrimType', type);
    yTest = mdl.predict(xTest);
  catch
    yTest = NaN(size(xTest, 1), 1);
  end
end
