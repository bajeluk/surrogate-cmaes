classdef OrigRatioUpdaterConstant < OrigRatioUpdater
  properties
    origParams 
    lastRatio
  end
  
  methods 
    % get new value of parameter
    function newRatio = update(obj, ~, ~, ~, ~, ~)
      % ratio is not updated
      newRatio = obj.lastRatio;
    end

    function obj = OrigRatioUpdaterConstant(parameters)
      % constructor
      obj = obj@OrigRatioUpdater(parameters);
      if isstruct(parameters)
        % starting value of ratio for initial generations
        obj.lastRatio = defopts(parameters, 'startRatio', 0.2);
      elseif iscell(parameters)
        parsedParams = struct(parameters{:});
        obj.lastRatio = defopts(parsedParams, 'startRatio', 0.2);
      else
        assert(isnumeric(parameters), 'Invalid parameter for OrigRatioUpdaterConstant');
        obj.lastRatio = parameters;
      end
    end
    
  end
end