classdef DTFileStatistics < Observer
%SCREENSTATISTICS -- print statistics from DoubleTrainEC on screen
  properties
    verbosity
    datapath
    exp_id
    instance
    expFileID
    file
  end

  methods
    function obj = DTFileStatistics(params)
      obj@Observer();
      obj.verbosity = defopts(params, 'verbose', 5);
      obj.datapath  = defopts(params, 'datapath', '/tmp');
      obj.exp_id    = defopts(params, 'exp_id', datestr(now,'yyyy-mm-dd_HHMMSS'));
      obj.instance  = defopts(params, 'instance', NaN);
      obj.expFileID = defopts(params, 'expFileID', '');
      obj.file  = [obj.datapath filesep obj.exp_id '_log_' obj.expFileID '.dat'];
    end

    function notify(obj, ec, varargin)
      % get the interesting data and process them
      fid = fopen(obj.file, 'a');
      if (ec.cmaesState.countiter == 1)
        % header lines
        % TODO: do not separate sections for independent restarts
        %       of the same instance
        fprintf(fid, '#  f/dim: %s  instance: %d  date: %s\n', strrep(regexprep('8_5D_5', '_[0-9]+$', ''), '_', ' '), ...
            obj.instance, datestr(now,'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, '#  iter  totEvals  orEvalsPop  Dopt  rmseReeval  rnkReeval  rnk2Models rnkValid  modelTypeUsed  nTrainData  nDataInRange\n');
      end
      model = 0;
      nTrainData = 0;
      if (~isempty(ec.model) && ec.model.isTrained() ...
          && ec.model.trainGeneration == ec.cmaesState.countiter)
        model = 1; nTrainData = ec.model.getTrainsetSize(); end
      if (~isempty(ec.retrainedModel) && ec.retrainedModel.isTrained() ...
          && ec.retrainedModel.trainGeneration == ec.cmaesState.countiter)
        model = 2; nTrainData = ec.retrainedModel.getTrainsetSize(); end

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
