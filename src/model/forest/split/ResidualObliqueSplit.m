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
      for i = 1:numel(discrimTypes)
        discrimType = discrimTypes{i};
        try
          model = fitcdiscr(obj.split_X, c, 'DiscrimType', discrimType);
        catch
          % singular covariance matrix
          pseudoDiscrimType = strcat('pseudo', discrimType);
          model = fitcdiscr(obj.split_X, c, 'DiscrimType', pseudoDiscrimType);
        end
        candidate = obj.splitCandidate;
        candidate.splitter = obj.createModelSplitter(model);
        candidate.gain = splitGain.get(candidate.splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
end