%LOADEFES -- loads EFE values from results *.mat files
%
% [efes, funs, dims, stats] = loadEFEs(path)
% Goes into the directory @path, loads corresponding 'scmaes_params.mat' file
% with definition of settings. Then loads each of the results *.mat file,
% takes the EFE values from the corresponding results and returns all such values
% in the @efes return paramter. Different statistics are returned in @stats return
% parameter.
function [efes, funs, dims, stats] = loadEFEs(path)
  TARGET = 1e-8;

  try
    cd(path);
    load('scmaes_params.mat');
  catch
    error('Error: scmaes_params.mat file not found');
  end

  % load basic parameters
  [~, ~, ~, nAlgs, totalIds] = getParamsFromIndex(1, bbParamDef, sgParamDef, cmParamDef);
  for i = 1:length(bbParamDef)
    switch (bbParamDef(i).name)
    case 'functions'
      funs = [bbParamDef(i).values{:}];
    case 'dimensions'
      dims = [bbParamDef(i).values{:}];
    % case 'instances'
    %   instances = [bbParamDef(i).values{:}];
    % case 'maxfunevals'
    %   maxfunevals_str = bbParamDef(i).values{1}
    end;
  end
  nInstances = 0;
  efes = cell(nAlgs,1);
  iAlg = 0;

  % for each ID
  for id = 1:totalIds
    % iAlg -- the id of algorithm settings
    iAlg = mod(iAlg, nAlgs) + 1;
    if (isempty(efes{iAlg}))
      efes{iAlg} = cell(length(funs), length(dims));
    end

    % load the parameters
    bbParams = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
    nInstances = max(nInstances, length(bbParams.instances));

    % these two for-cycles should run only once!
    for dim = bbParams.dimensions
      for ifun = bbParams.functions

        % take the indices of the dimension and function
        dimId = find(dims == dim);
        funId = find(funs == ifun);
        efes{iAlg}{funId, dimId} = [];

        % filename of the file with results
        expFileID = [num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
        resultsFile = [path filesep exp_id '_results_' expFileID '.mat'];
        if (~exist(resultsFile, 'file'))
          continue;
        end

        try
          % try to load results
          res = load(resultsFile);

          maxfunevals = eval(bbParams.maxfunevals);
          evals = res.exp_results.evals;
          fbests = res.exp_results.fbests;

          % calculate EFEs
          reached = (fbests < TARGET);
          efes{iAlg}{funId, dimId} = zeros(1, nInstances);
          efes{iAlg}{funId, dimId}(reached)  = evals(reached);
          efes{iAlg}{funId, dimId}(~reached) = maxfunevals * (1 + 1/9 * log( fbests(~reached) / TARGET ) );

          clear('res');

        catch
          warning(['Cannot load ' resultsFile '.']);
        end
      end
    end
  end

  % Calculate statistics

  stats = struct();
  stats.medians = zeros(length(funs)*length(dims), nAlgs);
  stats.uquart   = zeros(length(funs)*length(dims), nAlgs);
  stats.lquart   = zeros(length(funs)*length(dims), nAlgs);

  % calculate statistics for each algorithm settings
  for i = 1:nAlgs
    row = 0;
    for f = 1:length(funs)
      for d = 1:length(dims)
        row = row + 1;
        % medians and quartiles from 15 instances
        stats.medians(row, i) = median(efes{i}{f,d});
        stats.uquart(row, i)   = quantile(efes{i}{f,d}, 0.75);
        stats.lquart(row, i)   = quantile(efes{i}{f,d}, 0.25);
      end
    end
  end

  % ranking of the medians and quartiles among algoritm settings
  for r = 1:size(stats.medians, 1)
    stats.median_ranks(r, :) = ranking(stats.medians(r, :));
    stats.uquart_ranks(r, :)  = ranking(stats.uquart(r, :));
    stats.lquart_ranks(r, :)  = ranking(stats.lquart(r, :));
  end

  % summarize and sort the rankings
  stats.sum_median_ranks = sum(stats.median_ranks, 1);
  [stats.sort_median_ranks, stats.sort_median_idx] = sort(stats.sum_median_ranks);
  stats.sum_uquart_ranks = sum(stats.uquart_ranks, 1);
  [stats.sort_uquart_ranks, stats.sort_uquart_idx] = sort(stats.sum_uquart_ranks);
  stats.sum_lquart_ranks = sum(stats.lquart_ranks, 1);
  [stats.sort_lquart_ranks, stats.sort_lquart_idx] = sort(stats.sum_lquart_ranks);

  if (length(stats.sort_median_idx) >= 10)
    disp('The best 10 settings (ids upper, sum-ranks of f/d''s lower):');
    disp([stats.sort_median_idx(1:10); ...
          stats.median_ranks(:,stats.sort_median_idx(1:10))]);
    disp('');

    disp('ID''s of the best 10 settings:');
    stats.best_ids = [];
    for i = 1:size(stats.medians,1)
      stats.best_ids = [stats.best_ids (stats.sort_median_idx(1:10) + (i-1)*nAlgs)];
    end
    stats.best_ids = sort(stats.best_ids);
    disp(stats.best_ids);
  end
end
