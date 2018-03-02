classdef XGBoostModel < RandomForestModel
  
  properties %(Access = protected)
    rf_objectiveGrad1Func % 1st gradient of loss function
    rf_objectiveGrad2Func % 2nd gradient of loss function
  end
  
  methods
    function obj = XGBoostModel(modelOptions)
      % constructor
      modelOptions.boosting = true;
      modelOptions.treeFunc = defopts(modelOptions, 'treeFunc', ...
      @(x) GradientTreeModel(x));
      obj = obj@RandomForestModel(modelOptions);
      
      % model specific options
      obj.rf_shrinkage = defopts(modelOptions, 'rf_shrinkage', 0.1);
      obj.rf_objectiveGrad1Func = defopts(modelOptions, 'rf_objectiveGrad1Func', ...
        @(y, yPred) 2*(yPred-y));
      obj.rf_objectiveGrad2Func = defopts(modelOptions, 'rf_objectiveGrad2Func', ...
        @(y, yPred) repmat(2, size(y)));
    end
    
    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)      
      nFeatures = obj.rf_nFeaturesToSample;
      if nFeatures <= 0
        nFeatures = size(X, 2);
      end
      nRows = round(size(X, 1) * obj.rf_inBagFraction);
      if nRows <= 0
        nRows = size(X, 1);
      end
      
      yPred = zeros(size(y));
      obj.rf_trees = repmat(obj.treeTemplate, obj.rf_nTrees, 1);
      for iTree = 1:obj.rf_nTrees
        sample = struct;
        sample.features = datasample(1:size(X, 2), nFeatures, 2, 'Replace', false);
        sample.idx = datasample((1:size(X, 1))', nRows, 1, 'Replace', obj.rf_sampleWithReplacement);
        sample.X = X(sample.idx, sample.features);
        sample.y = y(sample.idx, :);
        obj.rf_trees(iTree).features = sample.features;
        obj.rf_trees(iTree).model = obj.rf_treeFunc(obj.rf_treeOptions);

        % fit using gradients
        sample.yPred = yPred(sample.idx, :);
        g = obj.rf_objectiveGrad1Func(sample.y, sample.yPred);
        h = obj.rf_objectiveGrad2Func(sample.y, sample.yPred);
        obj.rf_trees(iTree).model = obj.rf_trees(iTree).model.trainModel(...
            sample.X, [g h]);
        yPredNew = obj.rf_trees(iTree).model.modelPredict(X(:, sample.features));
        obj.rf_trees(iTree).weight = obj.rf_shrinkage;
        yPred = yPred + obj.rf_trees(iTree).weight * yPredNew;
      end
    end
  end
end