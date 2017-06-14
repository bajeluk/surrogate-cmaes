classdef (Abstract) OrigRatioUpdaterAbstractError < OrigRatioUpdater
% OrigRatioUpdaterAbstractError -- abstract class for adaptive origRatio updating
%   via model error, its exponential smoothing and linear transfer to origRatio
%
% Core method update() updates OrigRatioUpdater's state and returns the ratio
% of population-size which should be evaluated by the original fitness. The
% ratio is determined according to the previous results and actual model's
% error calculated by computeErr().
%
% Concrete subclasses implement different computeErr() methods which return
% (the best approximation of) current model error based on comparison between
% model- and original-evaluated points.
%
% This implementation uses
%
% (1) _exponential weighted average_ (or exponential smoothing in other words)
% of models' errors -- this smooth is specified by 'updateRate' parameter, and
% 'defaultErr' is used instead of smoothed value if no such is available
%
% (2) _linear transfer function_ from smoothed error value into ratio of point
% to be original evaluated -- this linear function is specified by parameters
% 'lowErr', 'highErr', 'minRatio' and 'maxRatio'.
%
% TODO
% [ ] make also transfer function to be replacable with something different
%
% =======
% = Log =
% =======
% 2016-10-19    - thorough documentation
%               - negative updates possibly different via 'updateRateDown' parameter

  properties
    surrogateOpts       % options, for DTS Updaters are DTAdaptive_*
    ec                  % reference to the corresponding EvolutionControl instance

    % actual state of Updater:
    lastRatio           % last calculated valid origRatio
    lastUpdateGeneration % generation number when last update was done
    gain                % last used 'gain' of linear transfer function

    % transfer function definition:
    lowErr              % error <= 'lowErr' is neglected, saturates to 'minRatio'
    highErr             % error >= 'highErr' saturates to 'maxRatio'
    minRatio            % origRatio set with zero or neglactable error
    maxRatio            % origRatio set with high error rates

    % exponential weighted average/smooth definition:
    updateRate          % update rate of the exponential smoothing (for positive updates)
    updateRateDown      % update rate for negative updates
                        % (current err is lower than last smoothed)
    defaultErr          % error-value used when no smoothed history value is found

    % history values of error and origRatio
    historyErr          % history values of calculated errors
    historySmoothedErr  % history values of EWA smoothed error values
    historyRatios       % history of used origRatios
  end

  methods (Abstract)
    err = computeErr(obj, modelY, origY, varargin);
  end

  methods
    % get new value of parameter
    function ratio = update(obj, modelY, origY, dim, lambda, countiter, varargin)
      % ratio is updated via
      % (1) exponential smooth (exp. weighted average) of model error 
      % (2) linear transfer from smoothed error into origRatio
      %
      % Notes:
      % - update() should be called every generation; if all the samples
      %   were evaluated with the original fitness, model should be also
      %   constructed and reasonable modelY values should be passed here
      % - if update() is not called in any particular generation(s),
      %   it results in NaN entries in historyErr for that generation(s)

      RATIO_TOLERANCE = 0.05;
      MAX_RATIO_ITERATES = 100;

      if (nargin >= 7), obj.ec = varargin{1}; end
      obj.historyErr((obj.lastUpdateGeneration+1):(countiter-1)) = NaN;
      obj.historySmoothedErr((obj.lastUpdateGeneration+1):(countiter-1)) = NaN;
      obj.historyRatios((obj.lastUpdateGeneration+1):(countiter-1)) = NaN;

      % compute current error and add it to history
      if (isempty(modelY) || (max(modelY) - min(modelY)) == 0 ...
          || (max(origY) - min(origY)) == 0)
        % there's not enough data to calculate current model's error
        err = NaN;
      else
        % calculate the error according to the subclass' computeErr()
        mu = ceil(obj.ec.cmaesState.mu * (size(modelY, 2) / obj.ec.cmaesState.lambda));
        err = obj.computeErr(modelY, origY, mu);
      end
      obj.historyErr(countiter) = err;

      % calculate exponentialy smoothed value of the error
      %      e_{t} = (1-a) * e_{t-1}  +  a * err
      % and add it to the history
      %
      lastIdx = find(~isnan(obj.historySmoothedErr(1:max(countiter-1,0))), 1, 'last');
      if (isempty(lastIdx))
        % there's no valid smoothed error value in history, use 'defaultErr'
        lastSmoothedErr = obj.defaultErr;
      else
        % take the last non-NaN smoothed error value
        lastSmoothedErr = obj.historySmoothedErr(lastIdx);
      end
      if (~isnan(err))
        % we have got a new current error ==> use EWA smooth
        if (err > lastSmoothedErr)
          % err is higher than last smoothed
          smoothedErr = (1-obj.updateRate) * lastSmoothedErr + obj.updateRate * err;
        else
          % err is lower than last smoothed
          smoothedErr = (1-obj.updateRateDown) * lastSmoothedErr + obj.updateRateDown * err;
        end

        obj.historySmoothedErr(countiter) = smoothedErr;
      else
        % we do not have a new error value, just use the last one
        smoothedErr = lastSmoothedErr;
        obj.historySmoothedErr(countiter) = NaN;
      end

      % transform the smoothed error into the ratio using linear transfer function
      % Note: the while cycle is due to the possible dependence of lowErr/highErr
      %       on the ratio, which should be the current ratio (not from the last generation)
      prevRatio = -1;
      ratio = obj.lastRatio;
      iterates = 0;
      while (abs(prevRatio - ratio) > RATIO_TOLERANCE)
        if (isa(obj.lowErr, 'function_handle'))
          lowErrValue = obj.lowErr([dim, ratio]);
        else
          lowErrValue = obj.lowErr;
        end
        if (isa(obj.highErr, 'function_handle'))
          highErrValue = obj.highErr([dim, ratio]);
        else
          highErrValue = obj.highErr;
        end
        obj.gain = min(max(0, smoothedErr - lowErrValue), ...
            (highErrValue - lowErrValue)) / (highErrValue - lowErrValue);
        prevRatio = ratio;
        ratio = obj.minRatio + obj.gain * (obj.maxRatio - obj.minRatio);

        iterates = iterates + 1;
        if (iterates >= MAX_RATIO_ITERATES)
          warning('OrigRatioUpdate: ratio calculation diverges!');
          break;
        end
      end

      % save the calculated ratio into history
      obj.lastRatio = ratio;
      obj.historyRatios(countiter) = ratio;
      obj.lastUpdateGeneration = countiter;
    end

    function obj = OrigRatioUpdaterAbstractError(ec, parameters)
      % constructor
      obj = obj@OrigRatioUpdater();
      obj.ec = ec;
      obj.surrogateOpts = parameters;

      % default parameter settings
      obj.lowErr     = myeval(defopts(obj.surrogateOpts, 'DTAdaptive_lowErr',  0.15));
      obj.highErr    = myeval(defopts(obj.surrogateOpts, 'DTAdaptive_highErr', 0.40));
      obj.updateRate = defopts(obj.surrogateOpts, 'DTAdaptive_updateRate', 0.40);
      obj.maxRatio   = defopts(obj.surrogateOpts, 'DTAdaptive_maxRatio',   1.00);
      obj.minRatio   = defopts(obj.surrogateOpts, 'DTAdaptive_minRatio',   0.05);
      obj.defaultErr = myeval(defopts(obj.surrogateOpts, 'DTAdaptive_defaultErr', 0.20));
      obj.lastRatio  = obj.minRatio;    % start conservatively...

      % for negative updates, use the same rate as positive rate as default;
      % this default is used also if [] is supplied in experiment definition
      % for 'updateRateDown'
      obj.updateRateDown = myeval(defopts(obj.surrogateOpts, 'DTAdaptive_updateRateDown', obj.updateRate));

      % history initialization
      obj.historyErr = [];
      obj.historySmoothedErr = [];
      obj.historyRatios = [];
      obj.lastUpdateGeneration = 0;
    end


    function value = getLastRatio(obj)
      % call "mock" update if more than one generation passed since last update()
      if (obj.ec.cmaesState.countiter > obj.lastUpdateGeneration + 1)
        obj.update([], [], [], [], obj.ec.cmaesState.countiter);
      end
      value = obj.lastRatio;
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
