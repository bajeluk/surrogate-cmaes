classdef DTICLogger < Observer
  %DTICLOGGER -- save log of information criteria for trained models
    
  properties
    datapath
    exp_id
    instance
    expFileID
    file
  end
  
  methods
    function obj = DTICLogger(params)
      obj@Observer()
      obj.datapath  = defopts(params, 'datapath', '/tmp');
      obj.exp_id    = defopts(params, 'exp_id', datestr(now,'yyyy-mm-dd_HHMMSS'));
      obj.instance  = defopts(params, 'instance', NaN);
      obj.expFileID = defopts(params, 'expFileID', '');
      obj.file  = [obj.datapath filesep obj.exp_id '_iclog_' obj.expFileID '.mat'];
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

        for f_cell = fieldnames(mdl.modelsIC)
          fn = f_cell{:};
          iclog.(fn) = [];
        end

        iclog.header = {'countiter', 'trainTrial', mdl.modelNames};
      end

      iclog.isTrained = [iclog.isTrained; ...
        [countiter obj.trainTrial mdl.isTrained]];

      for f_cell = fieldnames(mdl.modelsIC)
        fn = f_cell{:};
        iclog.(fn) = [iclog.(fn); ...
          [countiter obj.trainTrial mdl.modelsIC.(fn)]];
      end

      save(obj.file, '-struct', 'iclog');
  end
end

