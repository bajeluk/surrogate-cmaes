classdef CombiPolynomialModel < WeakModel
  
  properties %(Access = protected)
    weak_modelSpec % model specification from MATLAB fitlm function
    % (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    % except 'polyijk' settings
    weak_models   % individual models
    weak_useModel % id of model to use
    weak_coeff    % coefficients
    % weak_coeffCov % coefficient covariance
    % weak_features % used features
  end
  
  methods
    function obj = CombiPolynomialModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.weak_modelSpec = defopts(modelOptions, 'weak_modelSpec', 'constant');
      modelSpec_types = {'constant', 'linear', 'interactions', 'purequadratic', 'quadratic'};
      if ~iscell(obj.weak_modelSpec)
        obj.weak_modelSpec = {obj.weak_modelSpec};
      end
      for i = 1:numel(obj.weak_modelSpec)
        assert(any(strcmp(obj.weak_modelSpec{i}, modelSpec_types)), ...
          'Model ''%s'' cannot be specified as weak_modelSpec property for CombiPolynomialModel', ...
          obj.weak_modelSpec{i});
      end
      obj.weak_useModel = defopts(modelOptions, 'weak_useModel', 1:numel(obj.weak_modelSpec));
      obj.weak_coeff = defopts(modelOptions, 'weak_coeff', NaN);
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      [nPoints, dim] = size(X);
      for c = 1:numel(obj.weak_modelSpec)
        if ismember(c, obj.weak_useModel)
          % model is supposed to be trained
          if strcmp(obj.weak_modelSpec{c}, 'constant')
            weak_model = ConstantModel(struct('weak_coeff', obj.weak_coeff));
          else
            weak_model = RegressPolynomialModel(struct('weak_modelSpec', obj.weak_modelSpec{c}));
          end
          if weak_model.getMinTrainPoints(dim) > nPoints
            % not enough point to train the model
            obj.weak_models{c} = {};
          else
            obj.weak_models{c} = weak_model.trainModel(X, y);
          end
        else
          % model not is supposed to be trained
          obj.weak_models{c} = {};
        end
      end
    end
    
    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      nModels = numel(obj.weak_useModel);
      nPoints = size(X, 1);
      if nModels == 1
        if isempty(obj.weak_models{obj.weak_useModel})
          yPred = NaN(nPoints, 1);
          sd2 = NaN(nPoints, 1);
          ci = NaN(nPoints, 2);
        else
          [yPred, sd2, ci] = obj.weak_models{obj.weak_useModel}.modelPredict(X);
        end
      else
        yPred = cell(1, nModels);
        sd2 = cell(1, nModels);
        ci = cell(1, nModels);
        for c = 1:numel(obj.weak_modelSpec)
          if ~isempty(obj.weak_models{c})
            [yPred{c}, sd2{c}, ci{c}] = obj.weak_models{c}.modelPredict(X);
          end
        end
      end
    end
    
    function N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
      if any(strcmp(obj.weak_modelSpec, 'constant'))
          N = ones(size(dim));
      elseif any(strcmp(obj.weak_modelSpec, 'linear'))
          N = 1 + dim;
      elseif any(strcmp(obj.weak_modelSpec, 'interactions'))
          N = 1 + dim + dim.*(dim-1)/2;
      elseif any(strcmp(obj.weak_modelSpec, 'purequadratic'))
          N = 1 + 2*dim;
      elseif any(strcmp(obj.weak_modelSpec, 'quadratic'))
          N = 1 + 2*dim + dim.*(dim-1)/2;
      else
          N = [];
      end
    end
    
    function obj = setUseModel(obj, id)
    % sets models to be used from weak_modelSpec
      assert(all(ismember(id, 1:numel(obj.weak_modelSpec))), ...
             'Some model ID is out of possible values range 1:%d', ...
             numel(obj.weak_modelSpec))
      obj.weak_useModel = id;
    end
    
  end
  
end

