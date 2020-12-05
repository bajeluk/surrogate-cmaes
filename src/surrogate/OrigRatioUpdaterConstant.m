classdef OrigRatioUpdaterConstant < OrigRatioUpdater
% OrigRatioUpdaterConstant -- ratio updater class where the origRatio 
%                             remains constant

  properties
    origParams 
    
    % actual state of Updater:
    lastRatio           % last calculated valid origRatio
    lastUpdateGeneration % generation number when last update was done
    gain                % last used 'gain'
    
    % history values of error and origRatio
    historyErr          % history values of calculated errors
    historySmoothedErr  % history values of EWA smoothed error values
    historyRatios       % history of used origRatios
  end
  
  methods 
    % get new value of parameter
    function newRatio = update(obj, ~, ~, ~, ~, countiter, varargin)
      % ratio is not updated
      newRatio = obj.lastRatio;
      
      % save the calculated ratio into history
      obj.historyErr((obj.lastUpdateGeneration+1):countiter) = NaN;
      obj.historySmoothedErr((obj.lastUpdateGeneration+1):countiter) = NaN;
      obj.historyRatios((obj.lastUpdateGeneration+1):countiter) = obj.lastRatio;
      obj.lastUpdateGeneration = countiter;
    end

    function obj = OrigRatioUpdaterConstant(parameters)
      % constructor
      obj = obj@OrigRatioUpdater();
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
      
      obj.gain = 0;
      % history initialization
      obj.historyErr = [];
      obj.historySmoothedErr = [];
      obj.historyRatios = [];
      obj.lastUpdateGeneration = 0;
    end
    
  end
end