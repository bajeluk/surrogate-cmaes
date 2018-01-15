classdef LmfitPolynomialModel < WeakModel
  
  properties %(Access = protected)
    weak_modelSpec % model specification for MATLAB fitlm function
    % (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    weak_model     % trained model
  end
  
  methods
    function obj = LmfitPolynomialModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.weak_modelSpec = defopts(modelOptions, 'weak_modelSpec', 'constant');
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      if size(X, 1) == 1
        X = [X; X];
        y = [y; y];
      end
      warning('off', 'stats:LinearModel:RankDefDesignMat');
      obj.weak_model = fitlm(X, y, obj.weak_modelSpec);
      warning('on', 'stats:LinearModel:RankDefDesignMat');
    end
    
    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      [yPred, ci] = obj.weak_model.predict(X);
      if nargout >= 2
        sd2 = confidenceToVar(ci);
      end
    end
    
    function N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
      switch obj.weak_modelSpec
        case 'constant'
          N = ones(size(dim));
        case 'linear'
          N = 1 + dim;
        case 'interactions'
          N = 1 + dim + dim.*(dim-1)/2;
        case 'purequadratic'
          N = 1 + 2*dim;
        case 'quadratic'
          N = 1 + 2*dim + dim.*(dim-1)/2;
        otherwise
          warning(['Minimal number of training points is not properly ', ...
            'defined for such ''weak_modelSpec'' property of LmfitPolynomialModel. ', ... 
            'The returned minimal number is 1 + dim.'])
          N = 1 + dim;
      end
    end
    
  end
  
end

