classdef NoneEC < EvolutionControl
  properties 
    model
  end
  
  methods 
    function obj = NoneEC()
    % constructor
      obj.model = [];
    end
    
    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of double trained evolution control
    
      lambda = cmaesState.lambda;
      countiter = cmaesState.countiter;
      
      archive = archive.save(arxvalid', fitness_raw', countiter);

      surrogateStats = NaN(1, 2);
      origEvaled = false(1, lambda);
    end
    
  end
  
end
