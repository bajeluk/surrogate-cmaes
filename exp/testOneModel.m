function [stats, models, y_models, varargout] = testOneModel(modelType, modelOpts, ds, nSnapshots, opts)
% modelTest(modelType, modelOpts, ds) compute statistics on chosen model
%
% Input:
%   modelType   - type of tested model acc. to ModelFactory | string
%   modelOpts   - options of the tested model | struct
%   ds          - dataset for testing | struct
%   nSnapshots  - the number of snapshots recorded as a dataset per CMA-ES run;
%                 there must be the same number of datasets in 'ds.testSetX'
%   opts        - other options and settings:
%       .statistics      -- names of statistics to calculate; default: {'mse'}
%                           | cell-array of strings
%       .snapshotsToTest -- array of snapshots integer indices to test
%                           default: [1:nSnapshots] | array of integers
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
  opts.snapshotsToTest = defopts(opts, 'snapshotsToTest', 1:nSnapshots);
  opts.trySecondModel = defopts(opts, 'trySecondModel', false);
  opts.tolX = defopts(opts, 'tolX', 1e-12);

  % prepare output fields
  for st = 1:length(opts.statistics)
    stats.(opts.statistics{st}) = NaN(1, nSnapshots);
  end
  models = cell(1, nSnapshots);
  y_models = cell(1, nSnapshots);
  models2 = cell(1, nSnapshots);
  y_models2 = cell(1, nSnapshots);

  % cycle through all snapshots (usually 10)
  for i = opts.snapshotsToTest
    [lambda, dim] = size(ds.testSetX{i});

    % prepare archive for this snapshot one generation shorter
    % than the genertion where the snapshot was done
    thisArchive = ds.archive.duplicate();
    thisArchive = thisArchive.restrictToGenerations(1:(ds.generations(i)-1));

    % prepare the testing population
    thisPopulation = Population(lambda, dim);
    thisPopulation = thisPopulation.addPoints(ds.testSetX{i}', zeros(size(ds.testSetY{i}')), ...
        ds.testSetX{i}', NaN(size(ds.testSetX{i}')), 0);

    if (~isfield(ds, 'models') || length(ds.models) < i || isempty(ds.models{i}))
      % there is no valid model in the provided dataset
      m = ModelFactory.createModel(modelType, modelOpts, ds.means{i});

      % Default options for generating trainsets
      modelOpts.trainsetType = defopts(modelOpts, 'trainsetType', 'nearest');
      modelOpts.trainsetSizeMax = defopts(modelOpts, 'trainsetSizeMax', 15*dim);
      modelOpts.trainRange = defopts(modelOpts, 'trainRange', 1.0);

      % prepare the training set if needed
      if (~isfield(m.options, 'trainsetType') ...
          ||  strcmpi(modelOpts.trainsetType, 'parameters'))
        [X_train, y_train] = thisArchive.getTrainsetData(modelOpts.trainsetType, ...
            myeval(modelOpts.trainsetSizeMax), ...
            ds.cmaesStates{i}.xmean', modelOpts.trainRange, ...
            ds.cmaesStates{i}.sigma, ds.cmaesStates{i}.BD, thisPopulation);
      else
        X_train = [];
        y_train = [];
      end

      if (isa(m,'ModelPool'))
        for j = (m.historyLength+1):-1:1;
          if (isempty(ds.testSetX{i, j}))
            fprintf('Empty ds.testSetX{%d, %d}, not training ModelPool',i, j);
            continue;
          end
          thisPopulation = Population(lambda, dim);
          thisPopulation = thisPopulation.addPoints(ds.testSetX{i, j}', zeros(size(ds.testSetY{i, j}')), ...
          ds.testSetX{i,j}', NaN(size(ds.testSetX{i,j}')), 0);

          thisArchive = ds.archive.duplicate();
          thisArchive = thisArchive.restrictToGenerations(1:(ds.generations(i)-j));
          m = m.train(X_train, y_train, ds.cmaesStates{i,j}, ds.sampleOpts{i,j}, ...
              thisArchive, thisPopulation);
        end
      else
        m = m.train(X_train, y_train, ds.cmaesStates{i}, ds.sampleOpts{i}, ...
            thisArchive, thisPopulation);
      end

    else
      % there IS some model in the dataset, so use it
      m = ds.models{i};
    end


    if m.isTrained()
      y = ds.testSetY{i};
      y_models{i} = m.predict(ds.testSetX{i});

      % calculate and save statistics
      for st = 1:length(opts.statistics)
        fname = opts.statistics{st};
        stats.(fname)(i) = predictionStats(y, y_models{i}, opts.statistics{st});
      end
      fprintf('Model (gen. # %3d, %3d pts) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
        ds.generations(i), size(m.getDataset_y(), 1), stats.mse(i), stats.kendall(i), stats.rde(i));

      if (opts.trySecondModel)
        % Try also the second (retrained) model
        if (isfield(ds, 'models2') && length(ds.models2) >= i && ~isempty(ds.models2{i}))
          m2 = ds.models2{i};
        else
          m2 = m;
          origPointsIdx = (ds.archive.gens == ds.generations(i));
          m2.dataset = mergeNewPoints(m2.dataset, m2, ds.archive.X(origPointsIdx,:), ...
              ds.archive.y(origPointsIdx), opts.tolX);
          m2Archive = ds.archive.duplicate();
          m2Archive = m2Archive.restrictToGenerations(1:(ds.generations(i)));
          m2 = m2.train(m2.dataset.X, m2.dataset.y, ds.cmaesStates{i}, ds.sampleOpts{i}, ...
              m2Archive, thisPopulation);
        end
        if (m2.isTrained())
          y_models2{i} = m2.predict(ds.testSetX{i});
        else
          y_models2{i} = [];
        end
        models2{i} = m2;
      end
    else
      fprintf('Model (gen. # %3d) is not trained\n', ds.generations(i));
    end

    models{i} = m;
  end

  if (opts.trySecondModel)
    varargout = {models2, y_models2};
  end
end

function dataset = mergeNewPoints(dataset, model, newPointsX, newPointsY, tolX)
  XT = ( (model.trainSigma * model.trainBD) \ newPointsX')';

  for i = 1:size(newPointsX, 1)
    d = bsxfun(@minus, dataset.X, XT(i,:));
    if (~all(abs(d) <= tolX))
      dataset.X = [dataset.X; XT(i, :)];
      dataset.y = [dataset.y; newPointsY(i)];
    end
  end
end
