%GPRMODEL -- Gaussian Process model for S-CMA-ES using Matlab GPs
%
% original from:   https://github.com/repjak/surrogate-cmaes/blob/ce14b4ec2079e310c4c465f4473debf6757095c2/src/model/GprModel.m
% Author: repjak <j.repicky@gmail.com>
% Date:   Wed Aug 3 22:14:34 2016 +0200
%
%     Add the GprModel (fitrgp).
%
classdef GprModel < Model
  % Implements GP with RegressionGP from MATLAB
  properties    % inherited from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was trained
    trainSigma            % sigma of the generation when the model was trained
    trainBD               % BD of the generation when the model was trained
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling

    % GprModel specific properties
    stdY                  % standard deviation of Y in training set, for output normalization
    cov                   % a struct with cov fcn, hyperparameters vector and noise std
    options
    hyp
    nErrors
    fitErr
    gprMdl                % a RegressionGP object
    logModel              % display model object after each training
  end

  properties (Access = protected)
    % a list of known strings identifiers of covariance functions accepted
    % by fitrgp
    covFcnType = {'squaredexponential', 'matern32', 'matern52', 'ardsquaredexponential', ...
      'ardmatern32', 'ardmatern52'}
  end

  methods (Access = public)
    function obj = GprModel(modelOptions, xMean)
      assert(size(xMean,1) == 1, 'GpModel (constructor): xMean is not a row-vector.');

      % modelOpts structure
      if (isempty(modelOptions))
        obj.options = struct();
      else
        obj.options = modelOptions;
      end

      % computed settings
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;

      % Statistics and ML Toolbox check
      if (~license('checkout', 'statistics_toolbox'))
        warning('GprModel: Statistics and Machine Learning Toolbox license not available. Model cannot be used');
      end

      obj.cov = defopts(obj.options, 'cov',  struct('fcn', 'ardsquaredexponential'));
      if (~any(ismember(obj.cov.fcn, obj.covFcnType)))
        % a function handle
        if (~isfield(obj.cov, 'hyp'))
          error('Hyperparameters must be specified for custom covariance functions');
        end
        obj.cov.fcn = myeval(obj.cov.fcn);
      end

      obj.options.normalizeY = defopts(obj.options, 'normalizeY', true);
      obj.options.normalizeX = defopts(obj.options, 'normalizeX', true);
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.gprMdl = [];
      obj.logModel = 0;
      obj.nErrors = 0;
      obj.fitErr = 0;
    end

    function nData = getNTrainData(obj)
      nData = 3 * obj.dim;
    end

    function trained = isTrained(obj)
      trained = (strcmpi(class(obj.gprMdl), 'RegressionGP')) & ~obj.fitErr;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      obj.dataset.X = X;
      obj.dataset.y = y;

      % normalize y if specified
      if (obj.options.normalizeY)
        obj.shiftY = mean(y);
        obj.stdY = std(y);
        y = (y - obj.shiftY) / obj.stdY;
      else
        obj.shiftY = 0;
        obj.stdY = 1;
      end

      % used in some myevaled strings
      dim = obj.dim;
      obj.fitErr = 0;

      try
        if (isfield(obj.cov, 'hyp') && isfield(obj.cov, 'sigma'))
          % both sigma and hyperparameter specified
          hyp = myeval(obj.cov.hyp);
          sigma = myeval(obj.cov.sigma);

          obj.gprMdl = fitrgp(obj.getDataset_X(), y, ...
            'FitMethod', 'exact', ...
            'PredictMethod', 'exact', ...
            'Sigma', sigma, ...
            'KernelFunction', obj.cov.fcn, ...
            'KernelParameters', hyp, ...
            'Standardize', obj.options.normalizeX ...
          );
        elseif (isfield(obj.cov, 'hyp'))
          % compute default sigma inside fitrgp
          hyp = myeval(obj.cov.hyp);

          obj.gprMdl = fitrgp(obj.getDataset_X(), y, ...
            'FitMethod', 'exact', ...
            'PredictMethod', 'exact', ...
            'KernelFunction', obj.cov.fcn, ...
            'KernelParameters', hyp, ...
            'Standardize', obj.options.normalizeX ...
          );
        elseif (isfield(obj.cov, 'sigma'))
          % compute default hyperparameters inside fitrgp
          sigma = myeval(obj.cov.sigma);

          obj.gprMdl = fitrgp(obj.getDataset_X(), y, ...
            'FitMethod', 'exact', ...
            'PredictMethod', 'exact', ...
            'Sigma', sigma, ...
            'KernelFunction', obj.cov.fcn, ...
            'Standardize', obj.options.normalizeX ...
          );
        else
          % default sigma and default hyperparameters
          obj.gprMdl = fitrgp(obj.getDataset_X(), y, ...
            'FitMethod', 'exact', ...
            'PredictMethod', 'exact', ...
            'KernelFunction', obj.cov.fcn, ...
            'Standardize', obj.options.normalizeX ...
          );
        end
        
        obj.trainGeneration = generation;
      catch err
        disp(getReport(err));
        obj.nErrors = obj.nErrors + 1;
        obj.fitErr = 1;
        obj.trainGeneration = -1;
      end

      if (obj.logModel)
        disp('Model:');
        disp(obj.gprMdl);
      end
    end

    function [ypred, ysd] = modelPredict(obj, X)
      % make prediction
      % @ypred      -- predicted response
      % @ysd        -- predicted standard deviation
      if (strcmpi(class(obj.gprMdl), 'RegressionGP'))
        [ypred, ysd] = predict(obj.gprMdl, X);
        ypred = ypred * obj.stdY + obj.shiftY;
        ysd = ysd * obj.stdY;
      else
        ypred = [];
        ysd = [];
        warning('Model not trained');
      end
    end

    function gprMdl = getGpr(obj)
      % get gprMdl object
      gprMdl = obj.gprMdl;
    end
  end
end

function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
