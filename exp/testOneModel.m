function [stats, models, y_models] = testOneModel(modelType, modelOpts, ds, nSnapshots, opts)
% modelTest(modelType, modelOpts, ds) compute statistics on chosen model
%
% Input:
%   modelType   - type of tested model acc. to ModelFactory | string
%   modelOpts   - options of the tested model | struct
%   ds          - dataset for testing | struct
%   nSnapshots  - the number of snapshots recorded as a dataset per CMA-ES run
%   opts        - other options and settings (e.g.  opts.statistics)
%
% Output:
%   stats       - cell array of structures with calculated statistics (one per snapshot)
%   models      - trained models (one per snapshot) | cell-array
%   y_models    - predicted values (one per snapshot) | cell-array
%
% See Also:
%   datasetFromInstance
%
% E.g.
% 
%   [new_stats, models(i_data, :), y_models(i_data, :)] = ...
%       testOneModel(modelType{m}, modelOptions{m}, ...
%       data{f_data, d_data, i_data}, dataNSnapshots, opts);

  if (isempty(nSnapshots))
    nSnapshots = length(ds.testSetX);
  end
  opts.statistics = defopts(opts, 'statistics', { 'mse' });

  % prepare output fields
  for st = 1:length(opts.statistics)
    stats.(opts.statistics{st}) = NaN(1, nSnapshots);
  end
  models = cell(1, nSnapshots);
  y_models = cell(1, nSnapshots);

  % cycle through all snapshots (usually 10)
  for i = 1:nSnapshots

    m = ModelFactory.createModel(modelType, modelOpts, ds.means{i});
    m = m.train(ds.trainSetX{i}, ds.trainSetY{i}, ds.cmaesStates{i}, ds.sampleOpts);

    if m.isTrained()
      y = ds.testSetY{i};
      y_models{i} = m.predict(ds.testSetX{i});
      % calculate and save statistics
      stats.mse(i)      = predictionStats(y, y_models{i}, 'mse');
      stats.mzoe(i)     = predictionStats(y, y_models{i}, 'mzoe');
      stats.kendall(i)  = predictionStats(y, y_models{i}, 'kendall');
      stats.rankmse(i)  = predictionStats(y, y_models{i}, 'rankmse');
      stats.rankmzoe(i) = predictionStats(y, y_models{i}, 'rankmzoe');
      stats.rde(i)      = predictionStats(y, y_models{i}, 'rde');
    end
    
    models{i} = m;

    fprintf('Model (gen. # %3d) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
      ds.generations(i), stats.mse(i), stats.kendall(i), stats.rde(i));
  end
end
