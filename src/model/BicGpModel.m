classdef BicGpModel < Model
  %BICMODEL Model that selects from multiple covariance function models
  %   and selects the best according to the Bayesian information criterion.
  
  properties

    % BicGpModel specific fields
  end
  
  methods (Access = public)
    function obj = BicGpModel(modelOptions, xMean)
      % constructor
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
      obj.trainLikelihood = Inf;
      obj.useShift  = defopts(obj.options, 'useShift', false);

      % Optimization Toolbox check
      obj.options.trainAlgorithm = defopts(obj.options, 'trainAlgorithm', 'fmincon');
      if (strcmpi(obj.options.trainAlgorithm, 'fmincon') ...
          && ~license('checkout', 'optimization_toolbox'))
        warning('GpModel: Optimization Toolbox license not available. Switching to minimize().');
        obj.options.trainAlgorithm = 'minimize';
      end

      % GPLAB initialization
      params = resetparams;
      params = setfunctions(params, ...
        'v_covSum', 2, ...
        'v_covProd', 2, ...
        'v_covMask', 1, ...
        'v_covADD', obj.dim ...
      );
    params = setterminals(params, ...
      @() terminal

    end
  end
  
end

function [nlZ, dnlZ] = linear_gp(linear_hyp, s_hyp, inf, mean, cov, lik, x, y, linear_hyp_start, const_hyp_idx)
  % extend the vector of parameters by constant (i.e. not optimized) elements
  % taken from the vector of initial values
  linear_hyp_start(~const_hyp_idx) = linear_hyp;
  linear_hyp = linear_hyp_start;

  hyp = rewrap(s_hyp, linear_hyp');
  [nlZ, s_dnlZ] = gp(hyp, inf, mean, cov, lik, x, y);
  dnlZ = unwrap(s_dnlZ)';
  dnlZ = dnlZ(~const_hyp_idx);
end

function [c, ceq] = nonlincons(x)
  % checks if the values x(1:(end-4)) are within 2.5 off median
  MAX_DIFF = 2.5;
  ceq = [];
  assert(size(x,2) == 1, 'Argument for nonlincons is not a vector');
  c = zeros(size(x));
  % test only for covariance parameters
  % TODO: are there always 4 more parameters?!
  c = abs(x(1:end-4) - median(x(1:end-4))) - MAX_DIFF;
end
