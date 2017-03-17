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
    [lambda, dim] = size(ds.testSetX{i});

    % prepare archive for this snapshot one generation shorter
    % than the genertion where the snapshot was done
    thisArchive = ds.archive.duplicate();
    thisArchive = thisArchive.restrictToGenerations(1:(ds.generations(i)-1));

    % prepare the testing population
    thisPopulation = Population(lambda, dim);
    thisPopulation = thisPopulation.addPoints(ds.testSetX{i}', zeros(size(ds.testSetY{i}')), ...
        ds.testSetX{i}', NaN(size(ds.testSetX{i}')), 0);

    % create the model
    m = ModelFactory.createModel(modelType, modelOpts, ds.means{i});

    % Default options for generating trainsets
    modelOpts.trainsetType = defopts(modelOpts, 'trainsetType', 'nearest');
    modelOpts.trainsetSizeMax = defopts(modelOpts, 'trainsetSizeMax', 15*dim);
    modelOpts.trainRange = defopts(modelOpts, 'trainRange', 1.0);

    % prepare the training set if needed
    if (~isfield(m.options, 'trainsetType') ...
        ||  strcmpi(modelOpts.trainsetType, 'parameters'))
      [X_train, y_train] = thisArchive.getTrainsetData('nearest', ...
          myeval(modelOpts.trainsetSizeMax), ...
          ds.cmaesStates{i}.xmean', modelOpts.trainRange, ...
          ds.cmaesStates{i}.sigma, ds.cmaesStates{i}.BD, thisPopulation);
    else
      X_train = [];
      y_train = [];
    end

    m = m.train(X_train, y_train, ds.cmaesStates{i}, ds.sampleOpts{i}, ...
        thisArchive, thisPopulation);

    if m.isTrained()
      y = ds.testSetY{i};
      y_models{i} = m.predict(ds.testSetX{i});

      % calculate and save statistics
      for st = 1:length(opts.statistics)
        fname = opts.statistics{st};
        stats.(fname)(i) = predictionStats(y, y_models{i}, opts.statistics{st});
      end
    end

    models{i} = m;

    fprintf('Model (gen. # %3d, %3d pts) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
      ds.generations(i), size(m.dataset.y, 1), stats.mse(i), stats.kendall(i), stats.rde(i));
  end
end
