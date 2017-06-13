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

      % prepare the population for 1st model training, with no points orig-evaluated
      thisPopulation = Population(lambda, dim);
      thisPopulation = thisPopulation.addPoints(ds.testSetX{i}', [], ...
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
          % prepare the population, with no points orig-evaluated
          thisPopulation = Population(lambda, dim);
          thisPopulation = thisPopulation.addPoints(ds.testSetX{i, j}', [], ...
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
    end


    if m.isTrained()
      % prepare archive just as it was in the model's generation after
      % its training
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

      % get the (first) model prediction
      y = ds.testSetY{i};
      y_models{i} = m.predict(ds.testSetX{i});

      % calculate and save statistics
      for st = 1:length(opts.statistics)
        fname = opts.statistics{st};
        stats.(fname)(i) = predictionStats(y, y_models{i}, opts.statistics{st});
      end
      % calculate RDE statistic on independent populations
      stats.rdeValid(i) = validationRDE(m, ds.cmaesStates{i}, 10, opts.bbob_func);

      if (isfield(stats, 'mse') && isfield(stats, 'kendall') && isfield(stats, 'rde'))
        fprintf('Model (gen. # %3d, %3d pts) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
            g, size(m.getDataset_y(), 1), stats.mse(i), stats.kendall(i), stats.rde(i));
      end

      if (opts.trySecondModel)
        % Try also the second (retrained) model
        if (isfield(ds, 'models2') && length(ds.models2) >= i && ~isempty(ds.models2{i}))
          % clone the loaded model
          m2 = ModelFactory.createModel(modelType, ds.models2{i}.options, ds.means{i});
          m2 = m2.clone(ds.models2{i});
        else
          % clone the first model (M1)
          m2 = ModelFactory.createModel(modelType, m.options, m.trainMean);
          m2 = m2.clone(m);
        end
        retrainM2 = ~ m2.isTrained();

        % add such points to the re-trained models' dataset which were
        % really orig-evaluated in this generation
        origPointsIdx = (thisArchive.gens == g);
        loadedM2DatasetSize = length(m2.dataset.y);
        m2.dataset = mergeNewPoints(m2.dataset, m2, thisArchive.X(origPointsIdx,:), ...
            thisArchive.y(origPointsIdx), opts.tolX);
        if (loadedM2DatasetSize < length(m2.dataset.y))
          % retrain the M2 model if the dataset has just enlarged
          retrainM2 = true;
        end

        nRequiredOrigPoints = ceil(opts.testOrigRatio * lambda);
        nRequiredNewOrigPoints = max(0, nRequiredOrigPoints - sum(origPointsIdx));

        if (nRequiredNewOrigPoints >= 1)
          % the testing origRatio is (probably) different to the used in
          % the experiment, so add new points into the population
          % Note: don't forget to erase 'model2' from the dataset 'ds'
          %       if you want the second model retrain!
          nOrigPointsInPop = length(thisPopulation.getOriginalY());
          newX = ds.testSetX{i}(nOrigPointsInPop + (1:nRequiredNewOrigPoints), :);
          newY = ds.testSetY{i}(nOrigPointsInPop + (1:nRequiredNewOrigPoints));
          phase = 1;
          thisPopulation = thisPopulation.updateYValue(newX', newY', nRequiredNewOrigPoints, phase);
          m2.dataset = mergeNewPoints(m2.dataset, m2, newX, newY, opts.tolX);
          retrainM2 = true;
        end

        if (retrainM2)
          % train the retrained model
          m2 = m2.train([], [], ds.cmaesStates{i}, ds.sampleOpts{i}, ...
              thisArchive, thisPopulation);
        end

        if (m2.isTrained())
          y_models2{i} = m2.predict(ds.testSetX{i});
          stats.rde2models(i) = errRankMu(y_models{i}, y_models2{i}, m2.stateVariables.mu);
          stats.rde2(i) = errRankMu(y_models2{i}, y, m2.stateVariables.mu);
          % calculate RDE statistic on independent populations
          stats.rdeValid2(i) = validationRDE(m2, ds.cmaesStates{i}, 10, opts.bbob_func);
          y_m1_replace = y_models{i};
          y_m1_replace(thisPopulation.origEvaled) = y(thisPopulation.origEvaled);
          stats.rdeM1_M1WReplace(i) = errRankMu(y_models{i}, y_m1_replace, m.stateVariables.mu);
          y_m2_replace = y_models2{i};
          y_m2_replace(thisPopulation.origEvaled) = y(thisPopulation.origEvaled);
          stats.rdeM1_M2WReplace(i) = errRankMu(y_models{i}, y_m2_replace, m.stateVariables.mu);
          stats.rdeM2_M2WReplace(i) = errRankMu(y_models2{i}, y_m2_replace, m2.stateVariables.mu);
        else
          fprintf('Model2 (gen. # %3d) is not trained\n', g);
          y_models2{i} = [];
          stats.rde2models(i) = NaN;
          stats.rde2(i) = NaN;
          stats.rdeValid2(i) = NaN;
          stats.rdeM1_M1WReplace(i) = NaN;
          stats.rdeM1_M2WReplace(i) = NaN;
          stats.rdeM2_M2WReplace(i) = NaN;
        end
        models2{i} = m2;
      end
    else
      fprintf('Model (gen. # %3d) is not trained\n', g);
      stats.rdeValid(i) = NaN;
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
    if (~all(min(abs(d), [], 1) <= tolX))
      dataset.X = [dataset.X; XT(i, :)];
      dataset.y = [dataset.y; newPointsY(i)];
    end
  end
end

function rde = validationRDE(model, cmaesState, nRepeats, bbob_func_handle)
  % RDE on multiple validation test sets
  rdeValid = NaN(1, nRepeats);
  for rep = 1:nRepeats
    [~, xValidTest, ~] = sampleCmaesNoFitness(model.trainSigma, cmaesState.lambda, ...
        cmaesState, model.sampleOpts);

    try
      preciseOpts = model.options;
      if (~isfield(preciseOpts, 'bbob_func'))
        preciseOpts.bbob_func = bbob_func_handle;
      end
      preciseModel = ModelFactory.createModel('bbob', preciseOpts, model.trainMean);
      yTest = preciseModel.predict(xValidTest');
      yPredict = model.predict(xValidTest');
    catch err
      warning('BBOB precise model cannot be used: %s', err.message);
      rde = NaN;
      return;
    end
    % kendallValid(i) = corr(yPredict, yTest, 'type', 'Kendall');
    % rmseValid(i) = sqrt(sum((yPredict - yTest).^2))/length(yPredict);
    rdeValid(rep) = errRankMu(yPredict, yTest, cmaesState.mu);
  end
  rde = nanmean(rdeValid);
end

