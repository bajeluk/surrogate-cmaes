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
  opts.testOrigRatio = defopts(opts, 'testOrigRatio', 0);

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
    if i > length(ds.testSetX)
      break;
    end
    
    [lambda, dim] = size(ds.testSetX{i});
    g = ds.generations(i);

    if (~isfield(ds, 'models') || length(ds.models) < i || isempty(ds.models{i}))
      % there is no valid model in the provided dataset
      m = ModelFactory.createModel(modelType, modelOpts, ds.means{i});

      % prepare archive for this snapshot one generation shorter
      % than the genertion where the snapshot was done
      thisArchive = ds.archive.duplicate();
      thisArchive = thisArchive.restrictToGenerations(1:(g-1));

      % prepare the testing population
      thisPopulation = Population(lambda, dim);
      thisPopulation = thisPopulation.addPoints(ds.testSetX{i}', zeros(size(ds.testSetY{i}')), ...
          ds.testSetX{i}', NaN(size(ds.testSetX{i}')), 0);

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
          thisArchive = thisArchive.restrictToGenerations(1:(g-j));
          m = m.train(X_train, y_train, ds.cmaesStates{i,j}, ds.sampleOpts{i,j}, ...
              thisArchive, thisPopulation);
        end
      else
        m = m.train(X_train, y_train, ds.cmaesStates{i}, ds.sampleOpts{i}, ...
            thisArchive, thisPopulation);
      end

    else
      % there IS some model in the dataset, so use it
      m = ModelFactory.createModel(modelType, ds.models{i}.options, ds.means{i});
      m = m.clone(ds.models{i});
      
      % prepare archive just as it was for the model
      thisArchive = ds.archive.duplicate();
      thisArchive = thisArchive.restrictToGenerations(1:(g));

      % prepare the testing population with the orig-evaluated points from
      % this generation
      thisPopulation = Population(lambda, dim);
      thisX = thisArchive.X(thisArchive.gens == g, :);
      thisY = thisArchive.y(thisArchive.gens == g);
      nOrigPoints = length(thisY);
      thisPopulation = thisPopulation.addPoints(thisX', thisY', ...
          thisX', NaN(size(thisX')), length(thisY), 0);
      % fill the rest of the population with the points from 'testSetX'
      % supplied from dataset 'ds'
      thisX = ds.testSetX{i}(nOrigPoints+1:end,:);
      thisPopulation = thisPopulation.addPoints(thisX', NaN(1,lambda - nOrigPoints), ...
          thisX', NaN(size(thisX')), 0, 4);      
    end


    if m.isTrained()
      % get the (first) model prediction
      y = ds.testSetY{i};
      y_models{i} = m.predict(ds.testSetX{i});

      % calculate and save statistics
      for st = 1:length(opts.statistics)
        fname = opts.statistics{st};
        stats.(fname)(i) = predictionStats(y, y_models{i}, opts.statistics{st});
      end
      if (isfield(stats, 'mse') && isfield(stats, 'kendall') && isfield(stats, 'rde'))
        fprintf('Model (gen. # %3d, %3d pts) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
            g, size(m.getDataset_y(), 1), stats.mse(i), stats.kendall(i), stats.rde(i));
      end

      if (opts.trySecondModel)
        % Try also the second (retrained) model
        if (isfield(ds, 'models2') && length(ds.models2) >= i && ~isempty(ds.models2{i}))
          m2 = ModelFactory.createModel(modelType, ds.models2{i}.options, ds.means{i});
          m2 = m2.clone(ds.models2{i});
        else
          m2 = m;
          % add such points to the re-trained models' dataset which were
          % really orig-evaluated in this generation
          origPointsIdx = (ds.archive.gens == g);
          m2.dataset = mergeNewPoints(m2.dataset, m2, ds.archive.X(origPointsIdx,:), ...
              ds.archive.y(origPointsIdx), opts.tolX);
          
          if (opts.testOrigRatio > 0)
            % the testing origRatio is (probably) different to the used in
            % the experiment, so add new points into the population
            % Note: don't forget to erase 'model2' from the dataset 'ds'
            %       if you want the second model retrain!
            nOrigPoints = opts.testOrigRatio * lambda;
            nNewOrigPoints = max(0, nOrigPoints - sum(thisArchive.gens == g));

            if (nNewOrigPoints > 1)
              nUsedOrigPoints = length(thisPopulation.getOriginalY());
              newX = ds.testSetX{i}(nUsedOrigPoints + (1:nNewOrigPoints), :);
              newY = ds.testSetY{i}(nUsedOrigPoints + (1:nNewOrigPoints));
              phase = 1;
              thisPopulation = thisPopulation.updateYValue(newX', newY', nNewOrigPoints, phase);
            end
          end

          % train the retrained model
          m2 = m2.train([], [], ds.cmaesStates{i}, ds.sampleOpts{i}, ...
              thisArchive, thisPopulation);
        end
        if (m2.isTrained())
          y_models2{i} = m2.predict(ds.testSetX{i});
          stats.rnk2Models(i) = errRankMu(y_models{i}, y_models2{i}, m2.stateVariables.mu);
        else
          y_models2{i} = [];
          stats.rnk2Models(i) = NaN;
        end
        models2{i} = m2;
      end
    else
      fprintf('Model (gen. # %3d) is not trained\n', g);
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
