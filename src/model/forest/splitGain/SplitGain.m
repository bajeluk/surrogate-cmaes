classdef (Abstract) SplitGain
% SplitGain evaluates split functions used in decision trees using
% custom metric. The evaluation can be based on the data alone or on the
% output of a model. If the model is polynomial, optimized evaluation is
% used.
  
  properties %(Access = protected)
    X % input data
    y % output data
    minSize % min size of one side of the split
    modelFunc % function which creates new model
    probabilistic % whether to compute also probability and confidence interval
    degree % degree of polynomial features
    XP % generated polynomial features
    allEqual % whether all y values are equal
    current % struct holding data for current node
    polyMethod % method for computing polynomial model - 'regress' or ''
  end
  
  methods
    function obj = SplitGain(options)
    % constructor
    % evaluates functions based on constant model
    %   (y, yPred, variancePred)
    % if 'degree' or 'modelFunc' is specified, uses custom model
    %   (y, yPred, variancePred)
    % if also 'probabilistic' is specified
    %   (y, yPred, variancePred, confidenceIntervalPred, probabilityPred)
    % options
    %   'minSize'       - specifies the min size of either side
    %   'degree'        - uses polynomial model of given degree
    %   'modelFunc'     - uses custom model
    %   'probabilistic' - computes probability and confidence interval
      if nargin >= 1
        obj.minSize = defopts(options, 'minSize', 1);
        obj.degree = defopts(options, 'degree', []);
        obj.polyMethod = defopts(options, 'polyMethod', '');
        obj.modelFunc = defopts(options, 'modelFunc', []);
        obj.probabilistic = defopts(options, 'probabilistic', false);
      end
    end
    
    function obj = reset(obj, X, y)
    % resets the evaluator with new data (X, y)
      [n, ~] = size(X);
      obj.X = X;
      obj.y = y;
      obj.allEqual = size(obj.y, 1) == 0 || all(obj.y == obj.y(1, :));
      if obj.allEqual
        return
      end
      if ~isempty(obj.degree)
        obj.XP = generateFeatures(X, obj.degree, true);
        obj.minSize = max(obj.minSize, size(obj.XP, 2));
      end
      obj.current = obj.getData(true(n, 1));
    end
    
    function gain = get(obj, splitter)
    % evaluates splitter function
      if obj.allEqual
        gain = -inf;
        return;
      end
      idx = splitter(obj.X);
      if sum(idx) < obj.minSize || sum(~idx) < obj.minSize
        gain = -inf;
        return;
      end
      left = obj.getData(idx);
      right = obj.getData(~idx);
      n = size(obj.y, 1);
      nLeft = size(left.y, 1);
      nRight = size(right.y, 1);
      gain = obj.current.value ...
        - nLeft/n * left.value ...
        - nRight/n * right.value;
      % gain = obj.current.value - left.value - right.value;
    end
  end
  
  methods (Abstract, Access = protected)
    value = getValue(obj, data)
    % evaluates data using custom metric
  end
  
  methods (Access = protected)
    function data = getData(obj, idx)
    % returns portion of data given by idx and stores its value
      data = struct;
      data.idx = idx;
      data.y = obj.y(idx, :);
      if isempty(obj.modelFunc) && isempty(obj.degree)
        % constant model
        coeff = mean(data.y);
        coeffCov = var(data.y);
        data.yPred = repmat(coeff, size(data.y));
        data.sd2 = repmat(coeffCov, size(data.y));
      elseif isempty(obj.modelFunc) && stricmp(obj.polyMethod, 'regress')
        % polynomial model using regress
        XP = obj.XP(idx, :);
        [coeff, ci] = regress(data.y, XP);
        features = coeff ~= 0;
        coeff = coeff(features);
        XP = XP(:, features);
        data.yPred = XP * coeff;
        data.sd2 = 
      elseif isempty(obj.modelFunc)
        % polynomial model
        XP = obj.XP(idx, :);
        M = XP' * XP;
        % check rank deficiency
        r = rank(M);
        if r < size(M, 2)
          % remove dependent columns
          [~, features] = rref(M);
          XP = XP(:, features);
          M = M(features, features);
        end
        warning('off', 'MATLAB:rankDeficientMatrix');
        warning('off', 'MATLAB:singularMatrix');
        warning('off', 'MATLAB:nearlySingularMatrix');
        Mi = inv(M);
        %coeff = Mi * XP' * data.y;
        %coeff = M \ (XP' * data.y);
        coeff = XP \ data.y;
        warning('on', 'MATLAB:rankDeficientMatrix');
        warning('on', 'MATLAB:singularMatrix');
        warning('on', 'MATLAB:nearlySingularMatrix');
        data.yPred = XP * coeff;
        % var(b) = E(b^2) * (X'*X)^-1
        % add realmin to avoid 0 * inf = NaN
        coeffCov = (mean((data.y - data.yPred).^2) + realmin) * Mi;
        data.sd2 = sum(XP * coeffCov .* XP, 2);
      else
        % custom model
        X = obj.X(idx, :);
        model = obj.modelFunc(xMean);
        model = model.trainModel(X, data.y, xMean, 0);
        if obj.probabilistic
          [data.yPred, data.sd2, data.ci, data.p] = ...
            model.modelPredict(X);
        else
          [data.yPred, data.sd2] = ...
            model.modelPredict(X);
        end
      end
      data.value = obj.getValue(data);
    end
  end
end