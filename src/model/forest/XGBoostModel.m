classdef XGBoostModel < RandomForestModel
  
  properties %(Access = protected)
    objectiveGrad1Func % 1st gradient of loss function
    objectiveGrad2Func % 2nd gradient of loss function
  end
  
  methods
    function obj = XGBoostModel(modelOptions, xMean)
      % constructor
      modelOptions.boosting = true;
      modelOptions.treeFunc = defopts(modelOptions, 'treeFunc', ...
      @(xMean) GradientTreeModel(struct, xMean));
      obj = obj@RandomForestModel(modelOptions, xMean);
      
      % model specific options
      obj.objectiveGrad1Func = defopts(modelOptions, 'objectiveGrad1Func', ...
        @(y, yPred) 2*(yPred-y));
      obj.objectiveGrad2Func = defopts(modelOptions, 'objectiveGrad2Func', ...
        @(y, yPred) 2);
    end
    
    function obj = trainModel(obj, X, y, xMean, generation)
      % train the model based on the data (X,y)
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      
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
        sample.xMean = mean(sample.X);
        obj.trees(iTree).features = sample.features;
        obj.trees(iTree).model = obj.treeFunc(sample.xMean);

        % fit using gradients
        sample.yPred = yPred(sample.idx, :);
        g = obj.objectiveGrad1Func(sample.y, sample.yPred);
        h = obj.objectiveGrad2Func(sample.y, sample.yPred);
        obj.trees(iTree).model = obj.trees(iTree).model.trainModel(...
            sample.X, [g h], sample.xMean, generation);
        yPredNew = obj.trees(iTree).model.modelPredict(X(:, sample.features));
        obj.trees(iTree).weight = obj.shrinkage;
        yPred = yPred + obj.trees(iTree).weight * yPredNew;
      end
    end
  end
end