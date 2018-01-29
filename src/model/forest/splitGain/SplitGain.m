classdef (Abstract) SplitGain
% SplitGain evaluates split functions used in decision trees using
% custom metric. The evaluation can be based on the data alone or on the
% output of a model. If the model is polynomial, optimized evaluation is
% used.
  
  properties %(Access = protected)
    splitGain_X % input data
    splitGain_y % output data
    splitGain_minSize % min size of one side of the split
    splitGain_modelFunc % function which creates new model
    splitGain_modelOpts % options for newly created model (if specified in splitGain_modelFunc)
    splitGain_weightedGain % whether gain is weighted by number of examples
    splitGain_degree % degree of polynomial features
    splitGain_XP % generated polynomial features
    splitGain_allEqual % whether all y values are equal
    splitGain_current % struct holding data for current node
    splitGain_polyMethod % method for computing polynomial model - 'regress' or ''
    splitGain_computeSd2 = false % computes sd2
    splitGain_computeMse = false % computes mse
  end
  
  methods
    function obj = SplitGain(options)
    % constructor
    % evaluates functions based on constant model
    % if 'splitGain_degree' or 'splitGain_modelFunc' is specified, uses 
    % custom model options
    %   'splitGain_minSize'       - specifies the min size of either side
    %   'splitGain_degree'        - uses polynomial model of given degree
    %   'splitGain_polyMethod'    - method for computing polynomial model
    %     ''            - default model
    %     'regress'     - regress function
    %   'splitGain_modelFunc'     - uses custom model
    %   'splitGain_weightedGain'  - whether gain is weighted by number of 
    %                               examples
      if nargin > 0
        obj.splitGain_minSize = defopts(options, 'splitGain_minSize', 1);
        obj.splitGain_degree = defopts(options, 'splitGain_degree', []);
        obj.splitGain_polyMethod = defopts(options, 'splitGain_polyMethod', '');
        obj.splitGain_modelFunc = defopts(options, 'splitGain_modelFunc', []);
        obj.splitGain_modelOpts = defopts(options, 'splitGain_modelOpts', options);
        obj.splitGain_weightedGain = defopts(options, 'splitGain_weightedGain', true);
      end
    end
    
    function obj = reset(obj, X, y)
    % resets the evaluator with new data (X, y)
      [n, dim] = size(X);
      obj.splitGain_X = X;
      obj.splitGain_y = y;
      varY = var(y);
      obj.splitGain_allEqual = any(isnan(varY)) || all(varY < eps(max(varY)) * n);
      if obj.splitGain_allEqual
        return
      end
      if ~isempty(obj.splitGain_degree)
        obj.splitGain_XP = generateFeatures(X, obj.splitGain_degree, true);
        obj.splitGain_minSize = max(obj.splitGain_minSize(dim), size(obj.splitGain_XP, 2));
      end
      obj.splitGain_current = obj.getData(true(n, 1));
    end
    
    function gain = get(obj, splitter)
    % evaluates splitter function
      if obj.splitGain_allEqual
        gain = -inf;
        return;
      end
      idx = splitter(obj.splitGain_X) <= 0.5;
      minSize = obj.splitGain_minSize(size(obj.splitGain_X, 2));
      if sum(idx) < minSize || sum(~idx) < minSize
        gain = -inf;
        return;
      end
      left = obj.getData(idx);
      right = obj.getData(~idx);
      n = size(obj.splitGain_y, 1);
      nLeft = size(left.y, 1);
      nRight = size(right.y, 1);
      if obj.splitGain_weightedGain
        gain = obj.splitGain_current.value ...
          - nLeft/n * left.value ...
          - nRight/n * right.value;
      else
        gain = obj.splitGain_current.value - left.value - right.value;
      end
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
      data.y = obj.splitGain_y(idx, :);
      if size(data.y, 2) > 1
        % gradient model
      elseif isempty(obj.splitGain_modelFunc) && isempty(obj.splitGain_degree)
        % constant model
        mu = sum(data.y) / numel(data.y);
        data.yPred = repmat(mu, size(data.y));
        if obj.splitGain_computeMse || obj.splitGain_computeSd2
          r = data.y - mu;
          data.mse = r' * r / numel(r);
        end
        if obj.splitGain_computeSd2
          data.sd2 = repmat(data.mse, size(data.y));
        end
      elseif isempty(obj.splitGain_modelFunc) && strcmpi(obj.splitGain_polyMethod, 'regress')
        % polynomial model using regress
        XP = obj.splitGain_XP(idx, :);
        warning('off', 'stats:regress:RankDefDesignMat');
        [~, ~, r, rint] = regress(data.y, XP);
        warning('on', 'stats:regress:RankDefDesignMat');
        data.yPred = data.y - r;
        ci = bsxfun(@plus, data.yPred - r, rint);
        if obj.splitGain_computeMse
          data.mse = r' * r / numel(r);
        end
        if obj.splitGain_computeSd2
          data.sd2 = confidenceToVar(ci);
        end
      elseif isempty(obj.splitGain_modelFunc)
        % polynomial model
        XP = obj.splitGain_XP(idx, :);
        warning('off', 'MATLAB:rankDeficientMatrix');
        warning('off', 'MATLAB:singularMatrix');
        warning('off', 'MATLAB:nearlySingularMatrix');
        data.yPred = XP * (XP \ data.y);
        if obj.splitGain_computeMse || obj.splitGain_computeSd2
          % var(b) = E(b^2) * (X'*X)^-1
          r = data.y - data.yPred;
          data.mse = r' * r / numel(r);
        end
        if obj.splitGain_computeSd2
          data.sd2 = data.mse * sum(XP / (XP' * XP) .* XP, 2);
        end
        warning('on', 'MATLAB:rankDeficientMatrix');
        warning('on', 'MATLAB:singularMatrix');
        warning('on', 'MATLAB:nearlySingularMatrix');
      else
        % custom model
        X = obj.splitGain_X(idx, :);
        model = obj.splitGain_modelFunc(obj.splitGain_modelOpts);
        model = model.trainModel(X, data.y);
        if obj.splitGain_computeSd2
          [data.yPred, data.sd2] = model.modelPredict(X);
        else
          [data.yPred] = model.modelPredict(X);
        end
        if obj.splitGain_computeMse
          if iscell(data.yPred)
            % multiple models
            for m = 1:numel(data.yPred)
              r = data.y - data.yPred{m};
              data.mse(m) = r' * r / numel(r);
            end
          else
            % one model
            r = data.y - data.yPred;
            data.mse = r' * r / numel(r);
          end
        end
      end
      data.value = obj.getValue(data);
    end
  end
end