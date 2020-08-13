classdef RandomForestModel < WeakModel
  properties (Constant, Access = protected)
    treeTemplate = struct(... % template for trees
        'model', [], ... % model
        'features', [], ... % used features in the model
        'weight', 0); % weight of the model in the resulting tree
  end
  
  properties %(Access = protected)
    rf_boosting              % whether boosting is enabled
    rf_inBagFraction         % percentage of the data to use for training
    rf_lossFunc              % loss function
    rf_nFeaturesToSample     % number of variables
    rf_nTrees                % number of trees
    rf_oobError              % out of bag error
    rf_sampleWithReplacement % whether to use sampling with replacement
    rf_shrinkage             % shrinkage parameter
    rf_treeFunc              % function which creates a new tree
    rf_treeOptions           % tree settings and options
    rf_trees                 % trained trees
  end
  
  methods
    function obj = RandomForestModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      
      % model specific options
      obj.rf_treeOptions = defopts(modelOptions, 'rf_treeOptions', modelOptions);
      obj.rf_treeFunc = defopts(modelOptions, 'rf_treeFunc', @(x) TreeModel(x));
      obj.rf_nTrees = defopts(modelOptions, 'rf_nTrees', 100);
      obj.rf_nFeaturesToSample = defopts(modelOptions, 'rf_nFeaturesToSample', -1);
      obj.rf_sampleWithReplacement = defopts(modelOptions, 'rf_sampleWithReplacement', true);
      obj.rf_inBagFraction = defopts(modelOptions, 'rf_inBagFraction', 1);
      obj.rf_boosting = defopts(modelOptions, 'rf_boosting', false);
      obj.rf_shrinkage = defopts(modelOptions, 'rf_shrinkage', 0.5);
      obj.rf_lossFunc = defopts(modelOptions, 'rf_lossFunc', @mseLossFunc);
    end
    
    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      
      [nRows, nFeatures] = obj.trainDataProperties(X);
      
      yPred = zeros(size(y));
      obj.rf_trees = repmat(obj.treeTemplate, obj.rf_nTrees, 1);
      % tree training loop
      for iTree = 1:obj.rf_nTrees
        % sample data for tree training
        sample = obj.createSample(X, y, nFeatures, nRows);
        % create tree
        obj.rf_trees(iTree).features = sample.features;
        obj.rf_trees(iTree).model = obj.rf_treeFunc(obj.rf_treeOptions);
        
        % LS boosting
        if obj.rf_boosting
          % fit to residuals
          sample.yPred = yPred(sample.idx, :);
          r = sample.y - sample.yPred;
          obj.rf_trees(iTree).model = obj.rf_trees(iTree).model.trainModel(...
              sample.X, r);
          % find the best weight (simplified gradient of objective function)
          yPredNew = obj.rf_trees(iTree).model.modelPredict(X(:, sample.features));
          w = 1;
          objective = obj.rf_lossFunc(y, yPred + w * yPredNew);
          improved = true;
          while improved
            improved = false;
            eps = 0.01;
            for w1 = [w * (1 - eps), w * (1 + eps)]
              objective1 = obj.rf_lossFunc(y, yPred + w1 * yPredNew);
              if objective1 < objective
                w = w1;
                objective = objective1;
                improved = true;
                break;
              end
            end
          end
          obj.rf_trees(iTree).weight = w * obj.rf_shrinkage;
          yPred = yPred + obj.rf_trees(iTree).weight * yPredNew;
        
        % bagging
        else
          obj.rf_trees(iTree).model = obj.rf_trees(iTree).model.trainModel(...
              sample.X, sample.y);
          obj.rf_trees(iTree).weight = 1 / obj.rf_nTrees;
        end
      end
    end
    
    function [y, sd2] = modelPredict(obj, X)
      nX = size(X, 1);
      nTrees = size(obj.rf_trees, 1);
      y = zeros(nX, nTrees);
      sd2 = zeros(nX, nTrees);
      for iTree = 1:nTrees
        if obj.rf_trees(iTree).weight == 0
          % this tree doesn't contribute
          continue
        end
        XS = X(:, obj.rf_trees(iTree).features);
        [y(:, iTree), sd2(:, iTree)] = obj.rf_trees(iTree).model.modelPredict(XS);
      end
      weights = [obj.rf_trees(:).weight]';
      y = y * weights;
      % var(w1 * y1 + w2 * y2) = w1^2 * var(y1) + w2^2 * var(y2) + 2*w1*w2 * cov(y1,y2)
      % we are omitting the covariance between trees
      % TODO can we do this?
      sd2 = sd2 * (weights.^2);
    end
    
    function N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
      sampleTree = obj.rf_treeFunc(obj.rf_treeOptions);
      nFtToSample = obj.getNFeaturesToSample(dim);
      if nFtToSample > 0
        nFeatures = min(nFtToSample, dim);
      else
        nFeatures = dim;
      end
      N = sampleTree.getMinTrainPoints(nFeatures);
    end
    
  end
  
  methods(Access = protected)
    
    function [nData, nDim] = trainDataProperties(obj, X)
    % returns number of data and number of dimensions (features) for forest
    % training
    
      [N, dim] = size(X);
      % check number of features
      nDim = obj.getNFeaturesToSample(dim);
      if ~ismember(nDim, 1:dim)
        nDim = dim;
      end
      
      % check number of data
      nData = round(N * obj.rf_inBagFraction);
      if ~ismember(nData, 1:nData)
        nData = N;
      end
    end
    
    function sample = createSample(obj, X, y, nFeatures, nRows)
    % creates sample of training data
      sample.features = datasample(1:size(X, 2), nFeatures, 2, 'Replace', false);
      sample.idx = datasample((1:size(X, 1))', nRows, 1, 'Replace', obj.rf_sampleWithReplacement);
      sample.X = X(sample.idx, sample.features);
      sample.y = y(sample.idx, :);
    end
    
    function nFeatures = getNFeaturesToSample(obj, dim)
    % get n features to sample
      
      % number of features is gained through evaluation of expression in
      % rf_nFeaturesToSample property
      if ischar(obj.rf_nFeaturesToSample)
        try
          % the following expression can contain dim value e.g.
          % 'ceil(dim/3)'
          nFeatures = eval(obj.rf_nFeaturesToSample);
          assert(isnumeric(nFeatures), 'Result of eval function is not numerical.')
        catch err
          warning('rf_nFeaturesToSample could not be evaluated due to the following error: %s', err.message)
          nFeatures = dim;
        end
      % number of features is numerical
      elseif isnumeric(obj.rf_nFeaturesToSample)
        nFeatures = obj.rf_nFeaturesToSample;
      % wrong number of features format
      else
        warning(['Result of rf_nFeaturesToSample evaluation is not numeric. '...
                 'Setting number of features equal to dimension(%d)'], dim)
        nFeatures = dim;
      end
    end
    
  end
  
end