function ft = feature_ela_levelset(X, y, settings)
% ft = FEATURE_ELA_LEVELSET(X, y, settings) returns ELA levelset features
% for dataset [X, y] according to settings. 
% 
% The dataset [X, y] is split into two classes by a specific treshold
% calculated using different quantiles. Linear and quadratic discriminant 
% analysis are used to predict whether the objective values y fall below 
% or exceed the calculated threshold. The extracted features are based on 
% the distribution of the resulting cross-validated mean misclassification 
% errors of each classifier.
%
% settings:
%   methods   - discriminant analysis methods | default: {'lda', 'qda', 'mda'} 
%   nfolds    - number of folds in cross-validation | default: 10
%   quantiles - quantiles to calculate tresholds | default: [0.1, 0.25,
%               0.5]
% 
% Features:
%   mmce_[method]_[quantile]       - mean misclassification error of 
%                                    appropriate method using defined 
%                                    quantile
%   [method1]_[method2]_[quantile] - ratio between mean misclassification 
%                                    errors of two methods using defined 
%                                    quantile

  if nargin < 3
    if nargin < 2
      help feature_ela_levelset
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end
  
  emptyInput = false;
  % less than two points cannot be divided to groups (in cvpartition)
  if size(X, 1) < 2 || numel(y) < 2
    X = [NaN; NaN];
    y = X;
    emptyInput = true;
  end
  
  % initialize
  qnt = defopts(settings, 'quantiles', [0.1, 0.25, 0.5]);
  classMethods = defopts(settings, 'methods', {'lda', 'qda', 'mda'});
  nFolds = defopts(settings, 'nfolds', 10);
  
  nData = numel(y);
  nQuant = numel(qnt);
  nClassMet = numel(classMethods);
  
  % check number of folds
  if nData < nFolds
    nFolds = nData;
  end
  % initialize random number generator to gain the same CV partitions (and
  % identical feature values for identical data)
  rng_seed = rng;
  rng(nData, 'twister')
  % create instances for n-fold CV
  cp = cvpartition(nData, 'KFold', nFolds);
  % quantile tresholds
  y_tresh = quantile(y, qnt);
  mmce = NaN(nQuant, nClassMet);
  % quantile loop
  for q = 1:nQuant
    % class labels according to quantiles
    y_class = (y < y_tresh(q));
    for cm = 1:nClassMet
      switch classMethods{cm}
        case 'lda'
          classFcn = @(xTrain, yTrain, xTest) oneFold(xTrain, yTrain, xTest, 'linear', 'fitcdiscr');
        case 'qda'
          classFcn = @(xTrain, yTrain, xTest) oneFold(xTrain, yTrain, xTest, 'quadratic', 'fitcdiscr');
        case 'mda'
          classFcn = @(xTrain, yTrain, xTest) oneFold(xTrain, yTrain, xTest, '', 'mda');
        otherwise
          error('No such classification method as %s available for ELA levelset features implemented', ...
                classMethods{cm})
      end
      % CV
      mmce(q, cm) = crossval('mcr', X, y_class, 'predfun', classFcn, 'partition', cp);
      % mean misclassification error  features
      mmceName = sprintf('mmce_%s_%d', classMethods{cm}, 100*qnt(q));
      if emptyInput
        ft.(mmceName) = NaN;
      else
        ft.(mmceName) = mmce(q, cm);
      end
    end
    % mmce combination features
    combId = nchoosek(1:nClassMet, 2);
    for ci = 1:size(combId, 1)
      combiName = sprintf('%s_%s_%02d', ...
                          classMethods{combId(ci, 1)}, ...
                          classMethods{combId(ci, 2)}, ...
                          100*qnt(q));
      if emptyInput
        ft.(combiName) = NaN;
      else
        ft.(combiName) = mmce(q, combId(ci, 1)) / mmce(q, combId(ci, 2));
      end
    end
  end
  
  % return random number generator settings to original value
  rng(rng_seed)
  
end

function yTest = oneFold(xTrain, yTrain, xTest, type, discr_anal_fcn)
% one fold function
  try
    if strcmp(discr_anal_fcn, 'fitcdiscr')
      mdl = fitcdiscr(xTrain, yTrain, 'DiscrimType', type);
      yTest = mdl.predict(xTest);
    elseif strcmp(discr_anal_fcn, 'mda')
      J = numel(unique(yTrain));

      % FIXME: yTrain might not be consecutive integers
      yTrain = int8(yTrain);

      assert(all(sort(unique(yTrain))' == 1:max(unique(yTrain))));

      c = arrayfun(@(j) sum(yTrain == j), 1:J);
      % number of clusters per each group
      % must be less than the size of the smallest group - 1
      R = min(min(c) - 1, 3) * ones(J, 1);
      nIter = 5; % no. of iterations
      jProb = ones(1, J) / J; % uniform prior
      mdl = MDAModel(size(xTest, 2), J, R, nIter, jProb);
      mdl.fit(xTrain, yTrain);
      yTest = mdl.predict(xTest, 'map');
      yTest = nominal(yTest == 2);
    end
  catch err
    % uncomment for debugging:
    % rep = getReport(err);
    % fprintf('Failed with error: %s\n', rep);
    yTest = NaN(size(xTest, 1), 1);
  end
end