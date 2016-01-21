classdef NoneEC < EvolutionControl
  properties 
    model
  end
  
  methods 
    function obj = NoneEC()
    % constructor
      obj.model = [];
    end
    
    function [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of double trained evolution control
    
      lambda = cmaesState.lambda;
      countiter = cmaesState.countiter;
      
      if (strcmpi(func2str(surrogateOpts.sampleFcn), 'samplecmaes'))
        [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(cmaesState, sampleOpts, lambda, counteval, varargin{:});
      else
        error('surrogateManager: the only sampling method without model is "sampleCmaes()"');
      end
      
      archive = archive.save(arxvalid', fitness_raw', countiter);

      surrogateStats = NaN(1, 2);
    end
    
  end
  
end