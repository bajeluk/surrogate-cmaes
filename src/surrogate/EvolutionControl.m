classdef (Abstract) EvolutionControl
  properties (Abstract)
    model
  end
  
  methods (Abstract)
    % run one generation of evolution control
    [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats] = runGeneration(obj, cmaesState, surrogateOpts, archive, varargin)
  end
  
end