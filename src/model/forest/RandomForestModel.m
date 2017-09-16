classdef RandomForestModel < WeakModel
  properties (Constant, Access = private)
    treeTemplate = struct(... % template for trees
        'model', [], ... % model
        'features', [], ... % used features in the model
        'weight', 0); % weight of the model in the resulting tree
  end
  
  properties %(Access = protected)
    sampleWithReplacement % whether to use sampling with replacement
    inBagFraction % percentage of the data to use for training
    nFeaturesToSample % number of variables
    treeFunc % function which creates a new tree
    nTrees % number of trees
    trees % trained trees
    oobError % out of bag error
    boosting % whether boosting is enabled
    shrinkage % shrinkage parameter
    lossFunc % loss function
  end
  
  methods
    function obj = RandomForestModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      
      % model specific options
      obj.treeFunc = defopts(modelOptions, 'treeFunc', @() TreeModel(struct));
      obj.nTrees = defopts(modelOptions, 'nTrees', 1);
      obj.nFeaturesToSample = defopts(modelOptions, 'nFeaturesToSample', -1);
      obj.sampleWithReplacement = defopts(modelOptions, 'sampleWithReplacement', true);
      obj.inBagFraction = defopts(modelOptions, 'inBagFraction', 1);
      obj.boosting = defopts(modelOptions, 'boosting', false);
      obj.shrinkage = defopts(modelOptions, 'shrinkage', 0.1);
      obj.lossFunc = defopts(modelOptions, 'lossFunc', @immse);
    end
    
    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      
      nFeatures = obj.nFeaturesToSample;
      if nFeatures <= 0
        nFeatures = size(X, 2);
      end
      nRows = round(size(X, 1) * obj.inBagFraction);
      if nRows <= 0
        nRows = size(X, 1);
      end
      
      yPred = zeros(size(y));
      obj.trees = repmat(obj.treeTemplate, obj.nTrees, 1);
      for iTree = 1:obj.nTrees
        sample = struct;
        sample.features = datasample(1:size(X, 2), nFeatures, 2, 'Replace', false);
        sample.idx = datasample((1:size(X, 1))', nRows, 1, 'Replace', obj.sampleWithReplacement);
        sample.X = X(sample.idx, sample.features);
        sample.y = y(sample.idx, :);
        obj.trees(iTree).features = sample.features;
        obj.trees(iTree).model = obj.treeFunc();
        if obj.boosting
          % fit to residuals
          sample.yPred = yPred(sample.idx, :);
          r = sample.y - sample.yPred;
          obj.trees(iTree).model = obj.trees(iTree).model.trainModel(...
              sample.X, r, sample.xMean, generation);
          % find the best weight (simplified gradient of objective function)
          yPredNew = obj.trees(iTree).model.modelPredict(X(:, sample.features));
          w = 1;
          objective = obj.lossFunc(y, yPred + w * yPredNew);
          improved = true;
          while improved
            improved = false;
            eps = 0.01;
            for w1 = [w * (1 - eps), w * (1 + eps)]
              objective1 = obj.lossFunc(y, yPred + w1 * yPredNew);
              if objective1 < objective
                w = w1;
                objective = objective1;
                improved = true;
                break;
              end
            end
          end
          obj.trees(iTree).weight = w * obj.shrinkage;
          yPred = yPred + obj.trees(iTree).weight * yPredNew;
        else
          obj.trees(iTree).model = obj.trees(iTree).model.trainModel(...
              sample.X, sample.y, sample.xMean, generation);
          obj.trees(iTree).weight = 1 / obj.nTrees;
        end
      end
    end
    
    function [y, sd2] = modelPredict(obj, X)
      nX = size(X, 1);
      nTrees = size(obj.trees, 1);
      y = zeros(nX, nTrees);
      sd2 = zeros(nX, nTrees);
      for iTree = 1:nTrees
        if obj.trees(iTree).weight == 0
          % this tree doesn't contribute
          continue
        end
        XS = X(:, obj.trees(iTree).features);
        [y(:, iTree), sd2(:, iTree)] = obj.trees(iTree).model.modelPredict(XS);
      end
      weights = [obj.trees(:).weight]';
      y = y * weights;
      % var(w1 * y1 + w2 * y2) = w1^2 * var(y1) + w2^2 * var(y2) + 2*w1*w2 * cov(y1,y2)
      % we are omitting the covariance between trees
      % TODO can we do this?
      sd2 = sd2 * (weights.^2);
    end
  end
  
end