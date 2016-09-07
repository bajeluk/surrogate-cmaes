classdef DTFileStatistics < Observer
%SCREENSTATISTICS -- print statistics from DoubleTrainEC on screen
  properties
    verbosity
    file
  end

  methods
    function obj = DTScreenStatistics(params)
      obj@Observer();
      verbosity = defopts(params, 'verbose', 5);
      % TODO:
      % fileName = [surrogateOpts.experimentPath '_' expFileID '.dat'];
      file      = defopts(params, 'file',  ['/tmp/dt_file_log_' datestr(now,'yyyy-mm-dd_HHMMSS') '.dat']);
    end

    function notify(obj, ec, varargin)
      % get the interesting data and process them
      fid = fopen(file, 'a');
      if (ec.cmaesState.countiter == 1)
        fprintf(fid, '#  \n');
        fprintf(fid, '#  iter  totEvals  orEvalsPop  Dopt  rmseReeval  rnkReeval  rnk2Models rnkValid  modelTypeUsed  nTrainData  nDataInRange\n');
      end
      model = '0';
      nTrainData = 0;
      if (~isempty(ec.model) && ec.model.isTrained() ...
          && ec.model.trainGeneration == ec.cmaesState.countiter)
        model = '1'; nTrainData = ec.model.getTrainsetSize(); end
      if (~isempty(ec.retrainedModel) && ec.retrainedModel.isTrained() ...
          && ec.retrainedModel.trainGeneration == ec.cmaesState.countiter)
        model = '2'; nTrainData = ec.retrainedModel.getTrainsetSize(); end

      fprintf(fid, '%4d  %5d  %2d  %.4e  %.2e  %.2f  %.2f  %.2f  %d  %2d  %2d\n', ...
          ec.cmaesState.countiter, ec.counteval, sum(ec.pop.origEvaled), ...
          ec.stats.fmin - ec.surrogateOpts.fopt, ...
          ec.stats.rmseReeval, ...
          ec.stats.rankErrReeval, ...
          ec.stats.rankErr2Models, ...
          ec.stats.rankErrValid, ...
          model, nTrainData, ec.stats.nDataInRange ...
          );
      fclose(fid);
    end
  end
end
