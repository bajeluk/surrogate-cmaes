classdef ObserverFactory
  methods (Static)

    % create observers according to the settings in
    %   surrogateOpts.observers -- names of classes
    %   surrogateOpts.observersParams -- struct with optional settings
    % and register EvolutionControl 'ec' to these observers
    function [ec, observers] = createObservers(ec, surrogateOpts)
      observerNames  = defopts(surrogateOpts, 'observers', {});
      observerParams = defopts(surrogateOpts, 'observersParams', {});
      nObservers = length(observerNames);
      observers  = cell(1, nObservers);

      for i = 1:nObservers
        observerName = lower(observerNames{i});
        % parameters to the Observers' contructors can be a struct
        % or an empty array
        params = [];
        if (~isempty(observerParams))
          params = observerParams{i};
        end

        params.isAdaptive = isfield(surrogateOpts, 'updaterType') ...
            && ~strcmpi(surrogateOpts.updaterType, 'none');

        switch observerName
          case 'dtscreenstatistics'
            observers{i} = DTScreenStatistics(params);
            ec = observers{i}.registerObservable(ec);
          case 'dtfilestatistics'
            params.datapath  = surrogateOpts.datapath;
            params.exp_id    = surrogateOpts.exp_id;
            params.expFileID = surrogateOpts.expFileID;
            params.instance  = surrogateOpts.instance;
            observers{i} = DTFileStatistics(params);
            ec = observers{i}.registerObservable(ec);
          case 'dtmodelsaver'
            params.datapath  = surrogateOpts.datapath;
            params.exp_id    = surrogateOpts.exp_id;
            params.expFileID = surrogateOpts.expFileID;
            params.instance  = surrogateOpts.instance;
            observers{i} = DTModelSaver(params);
            ec = observers{i}.registerObservable(ec);
        end
      end

    end
  end
end

