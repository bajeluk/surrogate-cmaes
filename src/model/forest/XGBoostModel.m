classdef XGBoostModel < RandomForestModel
  
  properties %(Access = protected)
    rf_objectiveGrad1Func % 1st gradient of loss function
    rf_objectiveGrad2Func % 2nd gradient of loss function
  end
  
  methods
    function obj = XGBoostModel(modelOptions)
      % constructor
      modelOptions.rf_boosting = true;
      modelOptions.rf_treeFunc = defopts(modelOptions, 'treeFunc', ...
      @(x) GradientTreeModel(x));
      obj = obj@RandomForestModel(modelOptions);
      
      % model specific options
      obj.rf_shrinkage = defopts(modelOptions, 'rf_shrinkage', 0.1);
      if obj.rf_shrinkage > 1 || obj.rf_shrinkage < 0
        obj.rf_shrinkage = 0.1;
      end
      obj.rf_objectiveGrad1Func = defopts(modelOptions, 'rf_objectiveGrad1Func', ...
        @(y, yPred) 2*(yPred-y));
      obj.rf_objectiveGrad2Func = defopts(modelOptions, 'rf_objectiveGrad2Func', ...
        @(y, yPred) repmat(2, size(y)));
    end
    
    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)      
      
      [nRows, nFeatures] = obj.trainDataProperties(X);
      
      yPred = zeros(size(y));
      obj.rf_trees = repmat(obj.treeTemplate, obj.rf_nTrees, 1);
      for iTree = 1:obj.rf_nTrees
        % sample data for tree training
        sample = obj.createSample(X, y, nFeatures, nRows);
        % create tree
        obj.rf_trees(iTree).features = sample.features;
        obj.rf_trees(iTree).model = obj.rf_treeFunc(obj.rf_treeOptions);

        % fit using gradients
        sample.yPred = yPred(sample.idx, :);
        g = obj.rf_objectiveGrad1Func(sample.y, sample.yPred);
        h = obj.rf_objectiveGrad2Func(sample.y, sample.yPred);
        obj.rf_trees(iTree).model = obj.rf_trees(iTree).model.trainModel(...
            sample.X, [y g h]);
        yPredNew = obj.rf_trees(iTree).model.modelPredict(X(:, sample.features));
        % scale newly added weights by rf_shringkage after each step of 
        % tree boosting
        obj.rf_trees(iTree).weight = obj.rf_shrinkage;
        yPred = yPred + obj.rf_trees(iTree).weight * yPredNew;
      end
    end
  end
end