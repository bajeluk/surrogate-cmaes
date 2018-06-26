classdef DTICLogger < Observer
  %DTICLOGGER -- save log of information criteria for trained models
    
  properties
    datapath
    exp_id
    instance
    expFileID
    file
    printICs
  end

  methods
    function obj = DTICLogger(params)
      obj@Observer()
      obj.datapath  = defopts(params, 'datapath', '/tmp');
      obj.exp_id    = defopts(params, 'exp_id', datestr(now,'yyyy-mm-dd_HHMMSS'));
      obj.instance  = defopts(params, 'instance', NaN);
      obj.expFileID = defopts(params, 'expFileID', '');
      obj.file  = [obj.datapath filesep obj.exp_id '_iclog_' obj.expFileID '.mat'];
      obj.printICs = defopts(params, 'printICs', false);
    end

    function notify(obj, ec, varargin)
      if (isempty(ec.model) || ~isa(ec.model, 'ModelSelector')) && ...
        (isempty(ec.retrainedModel) || ~isa(ec.retrainedModel, 'ModelSelector'))
        return;
      end

      countiter = ec.cmaesState.countiter;

      if (~isempty(ec.retrainedModel))
        % contains both trials
        mdl = ec.retrainedModel;
      elseif (~isempty(ec.model))
        % the first model
        mdl = ec.model;
      end

      if (exist(obj.file, 'file'))
        iclog = load(obj.file);
      else
        iclog = struct();
        iclog.isTrained = [];
        iclog.nTrainErrors = [];
        iclog.bestIdx = [];
        iclog.mcmcConv = [];

        for f_cell = fieldnames(mdl.modelsIC)'
          fn = f_cell{:};
          iclog.(fn) = [];
        end

        iclog.header = {'countiter', 'trainTrial', mdl.modelNames};
      end

      if size(mdl.modelIsTrained,  1) == 1
        idx = [countiter mdl.trainTrial];
      else
        assert(mdl.trainTrial-1 == size(mdl.modelIsTrained, 1));
        idx = [repmat(countiter, mdl.trainTrial-1, 1) (1:(mdl.trainTrial-1))'];
      end

      iclog.isTrained = [iclog.isTrained; ...
        [idx mdl.modelIsTrained]];
      iclog.nTrainErrors = [iclog.nTrainErrors; ...
        [idx mdl.nTrainErrors]];
      iclog.bestIdx = [iclog.bestIdx; ...
        [idx mdl.bestIdx']];
      if isfield(mdl.modelsIC, 'rhat')
        converged = cellfun(@(r) all(r < 1.1), mdl.modelsIC.rhat, ...
            'UniformOutput', true);
        iclog.mcmcConv = [iclog.mcmcConv; ...
          [idx converged]];
      end

      for f_cell = fieldnames(mdl.modelsIC)'
        fn = f_cell{:};
        if ~iscell(mdl.modelsIC.(fn))
          iclog.(fn) = [iclog.(fn); ...
            [idx mdl.modelsIC.(fn)]];
        else
          iclog.(fn) = [iclog.(fn); ...
            [num2cell(idx) mdl.modelsIC.(fn)]];
        end
      end

      save(obj.file, '-struct', 'iclog');

      if obj.printICs
        ics = mdl.modelsIC.(mdl.ic);
        [~, idx] = sort(ics, 2, 'ascend');
        % calculate ranks
        ranks = zeros(size(ics));
        for i = 1:size(ics, 1)
          for j = 1:size(ics, 2)
            ranks(i, idx(i, j)) = j;
          end
        end

        for k = 1:size(ranks, 1)
          fprintf('[%6s] [%6s] %s\n', mdl.ic, mdl.modelNames{mdl.bestIdx(k)}, num2str(ranks(k, :)));
        end

        if isfield(mdl.modelsIC, 'rhat')
          rhatmax = cellfun(@(r) max(r), mdl.modelsIC.rhat, ...
            'UniformOutput', true);
          for k = 1:size(converged, 1)
            fprintf('[max(r)] %s\n', num2str(rhatmax(k, :)));
          end
        end
      end

    end
  end
end
