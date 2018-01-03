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
      % linear regression
      model = PolynomialModel(struct('modelSpec', obj.split_degree));
      model = model.trainModel(obj.split_X, obj.split_y);
      c = model.modelPredict(obj.split_X) < obj.split_y;
      
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
      best = obj.getDiscrAnal(splitGain, c, best, discrimTypes);
      
    end
  end
end