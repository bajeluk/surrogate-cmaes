function [stats, models, y_models, varargout] = testOneModel(modelType, modelOpts, ds, nSnapshots, opts)
% testOneModel(modelType, modelOpts, ds, nSnapshots, opts) compute 
% statistics on chosen model.
%
% [stats, models, y_models] = testOneModel(...) return statistics, trained
%   models and model predicted values.
%
% [..., models2, y_models2] = testOneModel(...) return also second training
%   models and predicted values.
%
% [..., models2, y_models2, modelOutput] = testOneModel(...) return also 
%   second training models, their predicted values, and outputs of first
%   training models.
%
% Input:
%   modelType   - type of tested model acc. to ModelFactory | string
%   modelOpts   - options of the tested model | struct
%   ds          - dataset for testing | struct
%   nSnapshots  - the number of snapshots recorded as a dataset per CMA-ES 
%                 run - there must be the same number of datasets in 
%                 'ds.testSetX' | integer
%   opts        - other options and settings:
%     .snapshotsToTest - array of snapshots integer indices to test;
%                         default: [1:nSnapshots] | array of integers
%     .statistics      - names of statistics to calculate; default: {'mse'}
%                        | cell-array of strings
%     .testModelOutput - test also output of the model (e.g. PoI, EI);
%                        default: false | boolean
%     .testOrigRatio   - value of origRatio (from DTS) used in tests;
%                        default: 0 | non-negative double
%     .tolX            - tolerance of distance between individual points;
%                        default: 1e-12 | non-negative double
%     .trySecondModel  - test also the second (retrained) model; default:
%                        false | boolean
%
% Output:
%   stats        - structures with calculated statistics (one per snapshot)
%                  | cell-array
%   models       - trained models (one per snapshot) | cell-array
%   y_models     - predicted values (one per snapshot) | cell-array
%   models2      - trained (second-training) models (one per snapshot) | 
%                  cell-array
%   y_models2    - predicted values of second-training models (one per 
%                  snapshot) | cell-array
%   modelOutputs - outputs of first-training models (one per snapshot) | 
%                  cell-array
%
% Example:
%
%   [new_stats, models(i_data, :), y_models(i_data, :)] = ...
%       testOneModel(modelType{m}, modelOptions{m}, ...
%       data{f_data, d_data, i_data}, dataNSnapshots, opts);
%
% See Also:
%   testModels, datasetFromInstances

  if (isempty(nSnapshots))
    nSnapshots = length(ds.testSetX);
  end
  opts.statistics = defopts(opts, 'statistics', { 'mse' });
  opts.snapshotsToTest = defopts(opts, 'snapshotsToTest', 1:nSnapshots);
  opts.trySecondModel = defopts(opts, 'trySecondModel', false);
  opts.tolX = defopts(opts, 'tolX', 1e-12);
  opts.testOrigRatio = defopts(opts, 'testOrigRatio', 0);
  opts.testModelOutput = defopts(opts, 'testModelOutput', false);

  % prepare output fields
  for st = 1:length(opts.statistics)
    stats.(opts.statistics{st}) = NaN(1, nSnapshots);
  end
  models = cell(1, nSnapshots);
  y_models = cell(1, nSnapshots);
  models2 = cell(1, nSnapshots);
  y_models2 = cell(1, nSnapshots);
  % prepare model outputs
  modelOutputs = cell(1, nSnapshots);
  if opts.testModelOutput
    if ~iscell(modelOpts.predictionType)
      modelOpts.predictionType = {modelOpts.predictionType};
    end
    for pt = 1:numel(modelOpts.predictionType)
      for st = 1:numel(opts.statistics)
        stats.outputs.(modelOpts.predictionType{pt}).(opts.statistics{st}) ...
          = num2cell(NaN(1, nSnapshots));
      end
    end
  end

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
            ds.cmaesStates{i}.sigma, ds.BDs{i}, thisPopulation);
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

          cmaesState = getCmaesState(ds, i, j);

          m = m.train(X_train, y_train, cmaesState, ds.sampleOpts{i,j}, ...
              thisArchive, thisPopulation);
        end
      else
        cmaesState = getCmaesState(ds, i);

        m = m.train(X_train, y_train, cmaesState, ds.sampleOpts{i}, ...
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
      cmaesState = getCmaesState(ds, i);
      stats.rdeValid(i) = validationRDE(m, cmaesState, 10, opts.bbob_func);

      if (isfield(stats, 'mse') && isfield(stats, 'kendall') && isfield(stats, 'rde'))
        fprintf('Model (gen. # %3d, %3d pts) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
            g, size(m.getDataset_y(), 1), stats.mse(i), stats.kendall(i), stats.rde(i));
      end
      
      % get model outputs
      if opts.testModelOutput
        % get outputs to cell-array
        modelOutputs{i} = m.getModelOutput(ds.testSetX{i});
        if ~iscell(modelOutputs{i})
          modelOutputs{i} = {modelOutputs{i}};
        end
        % model output type loop
        for mo = 1:numel(modelOutputs{i})
          if ~iscell(m.predictionType)
            predType = m.predictionType;
          else
            predType = m.predictionType{mo};
          end
          % sort output values
          if any(strcmpi(predType, {'sd2', 'poi', 'ei', 'erde', 'expectedrank'}))
            sortOrder = 'descend';
          else
            sortOrder = 'ascend';
          end
          [~, pointId] = sort(modelOutputs{i}{mo}, sortOrder);

          % calculate and save statistics for populations evaluated partially
          % by original fitness and partially by the model

          % loop according to the number of points evaluated using the 
          % original fitness
          for p = 1:numel(ds.testSetY{i})
            % create y-set composed of original and fitness evaluated
            % points
            ySet = y_models{i};
            ySet(pointId(1:p-1)) = y(pointId(1:p-1));
            % calculate statistics for this y-set
            for st = 1:length(opts.statistics)
              fname = opts.statistics{st};
              stats.outputs.(predType).(fname){i}(p) = predictionStats(y, ySet, fname);
            end
          end
        end
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
          cmaesState = getCmaesState(ds, i);
          m2 = m2.train([], [], cmaesState, ds.sampleOpts{i}, ...
              thisArchive, thisPopulation);
        end

        if (m2.isTrained())
          y_models2{i} = m2.predict(ds.testSetX{i});
          stats.rde2models(i) = errRankMu(y_models{i}, y_models2{i}, m2.stateVariables.mu);
          stats.rde2(i) = errRankMu(y_models2{i}, y, m2.stateVariables.mu);
          % calculate RDE statistic on independent populations
          stats.rdeValid2(i) = validationRDE(m2, getCmaesState(ds, i), 10, opts.bbob_func);
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

  if nargout > 3
    extraoutput = {models2, y_models2, modelOutputs};
    varargout = extraoutput(1:nargout-3);
  end
end

function dataset = mergeNewPoints(dataset, model, newPointsX, newPointsY, tolX)
  % Merge new points with the old points in the dataset taking into account
  % minimal coordinate difference 'tolX' (i.e., do not add identical, or
  % almost identical point to dataset).

  % no points to merge
  if isempty(newPointsX)
    return
  end
  % check numbers of new X and Y values
  assert(size(newPointsX, 1) == length(newPointsY), ...
         'scmaes:test1model:mergePointsUneq', ...
         'The number of X and Y values to merge with dataset does not match')

  % transform new points to sigma*BD space
  XT = ( (model.trainSigma * model.trainBD) \ newPointsX')';

  % no points in dataset -> add new points directly
  if isempty(dataset.X)
    dataset.X = XT;
    dataset.y = newPointsY;
  % some points in dataset -> add new points according to distances
  else
    for i = 1:size(newPointsX, 1)
      % calculate distances between new and dataset points in all
      % coordinates
      d = bsxfun(@minus, dataset.X, XT(i,:));
      % add new points distant further then tolerance
      if (~all(min(abs(d), [], 1) <= tolX))
        dataset.X = [dataset.X; XT(i, :)];
        dataset.y = [dataset.y; newPointsY(i)];
      end
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
    if isempty(yPredict) || any(isnan(yPredict))
      rdeValid(rep) = NaN;
    else
      rdeValid(rep) = errRankMu(yPredict, yTest, cmaesState.mu);
    end
  end
  rde = nanmean(rdeValid);
end


function cmaesState = getCmaesState(ds, varargin)
  % fill in fields that are stored in the outer level for memory
  % efficiency

  if nargin > 2
    i = varargin{1};
    j = varargin{2};
    cmaesState = ds.cmaesStates{i, j};
    cmaesState.BD = ds.BDs{i, j};
    cmaesState.diagD = ds.diagDs{i, j};
    cmaesState.diagC = ds.diagCs{i, j};
  elseif nargin > 1
    i = varargin{1};
    cmaesState = ds.cmaesStates{i};
    cmaesState.BD = ds.BDs{i};
    cmaesState.diagD = ds.diagDs{i};
    cmaesState.diagC = ds.diagCs{i};
  else
    error('Index not given');
  end
end
