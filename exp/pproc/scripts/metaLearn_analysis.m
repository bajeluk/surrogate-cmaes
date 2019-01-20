%% Read model options

exp_id = 'exp_DTSmodels_meta_02';
exppath_short = 'exp/experiments';
func = 1:24;
dim = [2 3 5 10 20];
inst = 11:15;
snModelIds = [1:6 8 9];
statistics = {'rde'};

exp_script_path = fullfile(exppath_short, [exp_id, '.m']);
run(exp_script_path);

opts.modelOptionsIndices = [1:6 8 9];

modelOptions_fullfact  = combineFieldValues(modelOptions);
modelOptionsTested = modelOptions_fullfact(opts.modelOptionsIndices);
nModels = length(modelOptionsTested);
modelOptionsTestedL = {'LIN', 'QUAD', 'SE', 'MAT5', 'RQ', 'NN', 'SE_PLUS_QUAD', 'GIBBS'};
assert(length(modelOptionsTestedL) == nModels, ...
  'Check whether model labels match model options.');

if exist(fullfile('exp/pproc', 'metaLearn_analysis_table.mat'), 'file')
  load(fullfile('exp/pproc', 'metaLearn_analysis_table.mat'));
else
  %% Read model results
  
  model_results = cell(1, nModels);
  
  missingData = zeros(0, 3);

  for m = 1:nModels
    modelHashName{m} = ['gpmodel_' modelHash(modelOptionsTested{m})];
    modelFolder{m} = fullfile(exppath_short, exp_id, [modelHashName{m}  '_250FE']);

    varNames = {'Function', 'Dimension', 'Instance', 'DataId', 'Snapshot', ...
      'mse', 'mzoe', 'kendall', 'rankmse', 'rankmzoe', 'rde', ...
      'rde2', 'rde2models', 'rdeValid', 'rdeValid2', 'rdeM1_M1WReplace', ...
      'rdeM1_M2WReplace', 'rdeM2_M2WReplace', 'mae', 'r2'};

    for c = 6:length(varNames)
      varNames{c} = [modelOptionsTestedL{m} '_' varNames{c}];
    end

    nCols = length(varNames);
    varTypes = repmat({'doublenan'}, 1, nCols);
    model_results{m} = table('Size', [0, nCols], 'VariableNames', varNames, 'VariableTypes', varTypes);

    for f = func
      for d = dim
        R = load(fullfile(modelFolder{m}, sprintf('%s_f%d_%dD.mat', modelHashName{m}, f, d)));

        for i = 1:size(R.stats.mse, 1)
          data_inst = find(R.instances == inst(i), 1);

          if isempty(data_inst)
            warning('f%d %dD inst%d missing\n', f, d, inst(i));
            missingData(end+1, :) = [opts.modelOptionsIndices(m), f, d];
            continue;
          end

          for modelId = 1:size(R.stats.mse, 2)
            data_modelId = find(snModelIds(modelId) == R.ids, 1);

            if isempty(data_modelId)
              warning('f%d %dD inst%d id%d missing\n', f, d, i, snModelIds(modelId));
              missingData(end+1, :) = [opts.modelOptionsIndices(m), f, d];
              continue;
            end

            for sn = 1:size(R.stats.mse, 3)
              model_results{m}(end+1, 1:nCols) = { ...
                f, d, data_inst, data_modelId, sn, ...
                R.stats.mse(i, modelId, sn), ...
                R.stats.mzoe(i, modelId, sn), ...
                R.stats.kendall(i, modelId, sn), ...
                R.stats.rankmse(i, modelId, sn), ...
                R.stats.rankmzoe(i, modelId, sn), ...
                R.stats.rde(i, modelId, sn), ...
                R.stats.rde2(i, modelId, sn), ...
                R.stats.rde2models(i, modelId, sn), ...
                R.stats.rdeValid(i, modelId, sn), ...
                R.stats.rdeValid2(i, modelId, sn), ...
                R.stats.rdeM1_M1WReplace(i, modelId, sn), ...
                R.stats.rdeM1_M2WReplace(i, modelId, sn), ...
                R.stats.rdeM2_M2WReplace(i, modelId, sn), ...
                R.stats.mae(i, modelId, sn), ...
                R.stats.r2(i, modelId, sn)};
            end
          end
        end
      end
    end
  end

  %% Read metafeatures
  mf_init = false;
  
  for f = func
    for d = dim
      fprintf('Loading MFs f%d %dD\n', f, d);
      for i = inst
        for modelId = snModelIds
          fname = fullfile(exppath_short, exp_id, 'metafeatures', ...
            sprintf('data_f%d_%dD_inst%d_id%d_fts.mat', f, d, i, modelId));

          if ~exist(fname, 'file'), continue; end
          Mfs=load(fname);

          if ~mf_init
            % parse metafeature names and initialize the table
            mf_names = strsplit(printStructure(Mfs.res(1).ft(1)));
            mf_names = mf_names(find(cellfun(@(s) contains(s, 'structure.'), mf_names)));
            mf_names = cellfun(@(x) ['mf_', replace(replace(x, 'structure.', ''), '.', '_')], mf_names, ...
              'UniformOutput', false);

            varNames = [{'Function', 'Dimension', 'Instance', 'DataId', 'Snapshot'}, ...
              mf_names];
            nCols = length(varNames);
            varTypes = repmat({'doublenan'}, 1, nCols);
            mf_results = table('Size', [0, nCols], 'VariableTypes', varTypes, 'VariableNames', varNames);
            mf_init = true;
          end

          for sn = 1:size(Mfs.res, 2)
            mf_results(end+1, 1:nCols) = [{f, d, i, modelId, sn}, ...
              num2cell(Mfs.res.values(:, sn)')];
          end
        end
      end
    end
  end

  results = mf_results;

  for m = 1:nModels
    results = outerjoin(results, model_results{m}, ...
      'Keys', {'Function', 'Dimension', 'Instance', 'DataId', 'Snapshot'}, ...
      'MergeKeys', true);
  end
  save(fullfile('/tmp', 'metaLearn_analysis_table'), 'model_results', 'mf_results', 'results');
end

%% Count NaNs for each model
stats_name = 'rde';
nan_counts = zeros(1, nModels);
for m = 1:nModels
  colName = [modelOptionsTestedL{m} '_' stats_name];
  nan_counts(m) = sum(isnan(model_results{m}.(colName)))/length(model_results{m}.(colName));
end

%% Friedman test on RDE
stats_name = 'rde';
colNames = cellfun(@(m) strjoin({m, stats_name}, '_'), modelOptionsTestedL, 'UniformOutput', false);
d = results{:, colNames};
% remove datasets where no models was trained
d(all(isnan(d), 2), :) = [];
% replace non-trained or no-prediction cases with infinite negative likelihood
d(isnan(d)) = Inf;
% run friedman
friedman(d)