classdef NoneEC < EvolutionControl & Observable
  properties
    model
    counteval
    surrogateOpts
    archive
    cmaesState
    stats
    pop
  end

  methods
    function obj = NoneEC(surrogateOpts)
    % constructor
      obj@Observable();
      obj.pop = [];
      obj.archive = [];
      obj.stats = struct();
      obj.model = [];
      obj.counteval = 0;
      obj.surrogateOpts = surrogateOpts;
      obj.stats = struct( ...
          'fmin', NaN ...    % minimal original fitness in population
      );
    end

    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of no evolution control

      lambda = cmaesState.lambda;
      dim = cmaesState.dim;
      countiter = cmaesState.countiter;
      obj.pop = Population(lambda, dim);
      obj.archive = archive;
      obj.cmaesState = cmaesState;

      [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(cmaesState, ...
          sampleOpts, lambda, counteval, 'Archive', archive, varargin{:});
      surrogateStats = NaN(1, 7);
      origEvaled = true(1, lambda);
      nonInfIdx = fitness_raw < Inf;
      obj.archive = archive.save(arxvalid(:,nonInfIdx)', fitness_raw(nonInfIdx)', countiter);
      obj.counteval = counteval;
      obj.pop = obj.pop.addPoints(arxvalid, fitness_raw, arx, arz, lambda, 2);

      % calculate statistics
      obj.stats.fmin = min(obj.pop.y);

      obj.notify_observers();
    end

  end

end
