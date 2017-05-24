classdef NoneEC < EvolutionControl & Observable
  properties 
    model
    counteval
  end

  methods 
    function obj = NoneEC()
    % constructor
      obj.model = [];
      obj.counteval = 0;
    end

    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of no evolution control

      lambda = cmaesState.lambda;
      countiter = cmaesState.countiter;

      [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(cmaesState, ...
          sampleOpts, lambda, counteval, 'Archive', archive, varargin{:});
      surrogateStats = NaN(1, 7);
      origEvaled = true(1, lambda);
      nonInfIdx = fitness_raw < Inf;
      archive = archive.save(arxvalid(:,nonInfIdx)', fitness_raw(nonInfIdx)', countiter);
      obj.counteval = counteval;
    end

  end

end
