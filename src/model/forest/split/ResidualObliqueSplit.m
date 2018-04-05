classdef ResidualObliqueSplit < Split
% ResidualObliqueSplit fits a polynomial model and splits the points into
% two classes based on which side of the line they lie

  properties %(Access = protected)
    split_degree % degree for fitted polynomial
  end
    
  methods
    function obj = ResidualObliqueSplit(options)
      obj = obj@Split(options);
      obj.split_degree = defopts(options, 'split_degree', 'constant');
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      % polynomial regression
      model = PolynomialModel(struct('weak_modelSpec', obj.split_degree));
      model = model.trainModel(obj.split_X, obj.split_y(:, 1));
      c = model.modelPredict(obj.split_X) < obj.split_y(:, 1);
      
      % discriminant analysis
      switch obj.split_degree
        case 'constant'
          discrimTypes = {'linear', 'quadratic'};
        case 'linear'
          discrimTypes = {'linear', 'quadratic'};
        case 'quadratic'
          discrimTypes = {'quadratic'};
        otherwise
          discrimTypes = {};
      end
      nDiscrTypes = numel(discrimTypes);
      [n, d] = size(obj.split_X);
      % get maximal number of hyperplanes
      maxHyp = obj.getMaxHyp(n, d);
      % lack of hyperplanes solve through random sampling from discriminant
      % types
      if maxHyp < nDiscrTypes
        discrimTypes = discrimTypes(randperm(nDiscrTypes, maxHyp));
      end
      best = obj.getDiscrAnal(splitGain, c, best, discrimTypes);
      
    end
  end
end