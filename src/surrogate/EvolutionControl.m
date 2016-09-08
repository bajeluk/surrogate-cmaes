classdef (Abstract) EvolutionControl
  % (handle object does not need to be returned in its own functions)
  properties (Abstract)
    model
    counteval
  end
  
  methods (Abstract)
    % run one generation of evolution control
    [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, archive, varargin)
  end
  
end
