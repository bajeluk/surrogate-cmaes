classdef DTModelSaver < Observer
%DTMODELSAVER -- save models from DoubleTrainEC into the *.mat file
  properties
    datapath
    exp_id
    instance
    expFileID
    file
  end

  methods
    function obj = DTModelSaver(params)
      obj@Observer();
      obj.datapath  = defopts(params, 'datapath', '/tmp');
      obj.exp_id    = defopts(params, 'exp_id', datestr(now,'yyyy-mm-dd_HHMMSS'));
      obj.instance  = defopts(params, 'instance', NaN);
      obj.expFileID = defopts(params, 'expFileID', '');
      obj.file  = [obj.datapath filesep obj.exp_id '_modellog_' obj.expFileID '.mat'];
    end

    % get the interesting data from DoubleTraineEC and process them
    function notify(obj, ec, varargin)
      if (exist(obj.file, 'file'))
        modellog = load(obj.file);
      else
        modellog = struct();
        modellog.models = cell(0);
        modellog.models2 = cell(0);
      end

      if (ec.cmaesState.countiter == 1 && ec.counteval <= (ec.cmaesState.lambda+1))
        modellog = struct();
        modellog.models = cell(0);
        modellog.models2 = cell(0);
        countiter = ec.cmaesState.countiter;
      elseif (exist(obj.file, 'file'))
        modellog = load(obj.file);
        countiter = length(modellog.models) + 1;
      else
        warning(['The mat-file ' obj.file ' with models should exist!']);
        countiter = ec.cmaesState.countiter;
      end

      if (~isempty(ec.model) && ec.model.isTrained() ...
          && ec.model.trainGeneration == ec.cmaesState.countiter)
        % the first model was successfully trained
        modellog.models{countiter} = ec.model;
      end
      if (~isempty(ec.retrainedModel) && ec.retrainedModel.isTrained() ...
          && ec.retrainedModel.trainGeneration == ec.cmaesState.countiter)
        % the second model was successfully trained
        modellog.models2{countiter} = ec.retrainedModel;
      end

      save(obj.file, '-struct', 'modellog');
    end
  end
end
