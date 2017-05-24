classdef NoneFileStatistics < Observer
%SCREENSTATISTICS -- log statistics into the file (for almost any EC)
  properties
    verbosity
    datapath
    exp_id
    instance
    expFileID
    file
    
    isHeaderWritten
  end

  methods
    function obj = NoneFileStatistics(params)
      obj@Observer();
      obj.datapath = defopts(params, 'datapath', '/tmp');
      obj.exp_id    = defopts(params, 'exp_id', datestr(now,'yyyy-mm-dd_HHMMSS'));
      obj.expFileID = defopts(params, 'expFileID', '');
      obj.file  = [obj.datapath filesep obj.exp_id '_log_' obj.expFileID '.dat'];
      obj.isHeaderWritten = false;
    end

    % get the interesting data from DoubleTraineEC and process them
    function notify(obj, ec, varargin)
      fid = fopen(obj.file, 'a');
      % write header
      if (ec.cmaesState.countiter == 1 && ec.counteval <= (ec.cmaesState.lambda+1))
        % only in the first iteration when counteval is very low
        % (do not write header after CMA-ES restarts)
        func_dim_id = strsplit(obj.expFileID, '_');
        fun = 0; dim = ec.cmaesState.dim; inst = func_dim_id{end};
        fprintf(fid, '#  f%d  dim: %dD  instance: %s  date: %s\n', ...
            fun, dim, inst, ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        colLabels = sprintf(['#  iter  totEvals  orEvalsPop  Fopt  sigmaPow2  ' ...
            'diagD_max  diagD_min  diagD_ratio  %s\n'], ...
            sprintf('diagD_%02d ', 1:length(ec.cmaesState.diagD)));
        fprintf(fid, colLabels);
        obj.isHeaderWritten = true;
      end

      %             iter  tot orEvl Dopt  sgm^2
      fprintf(fid, '%4d  %5d  %2d  %.4e  %.2e  ', ...
          ec.cmaesState.countiter, ec.counteval, sum(ec.pop.origEvaled), ...
          ec.stats.fmin, ec.cmaesState.sigma^2);
      fprintf(fid, ...
      ...% mxDi  miDi  Drio  D_1,...,D_n
          '%.2e  %.2e  %.2e  %s\n', ...
          max(ec.cmaesState.diagD), ...
          min(ec.cmaesState.diagD), ...
          max(ec.cmaesState.diagD)/min(ec.cmaesState.diagD), ...
          sprintf('%.2e ', sort(ec.cmaesState.diagD, 1, 'descend')) ...
          );
      fclose(fid);
    end
  end
end
