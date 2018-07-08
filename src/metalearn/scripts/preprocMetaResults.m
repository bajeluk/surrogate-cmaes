exppath = 'exp/experiments/exp_metaLearn_05';
S = load(fullfile(exppath, 'metalearn_params.mat'));
res_path = fullfile(exppath, 'results');
nModelTypes = length(S.modelParamDef);
res_fname_template = strjoin({'res_', ...
    '%dD_', ...
    'f%d_', ...
    'inst%d_', ...
    'N%d_', ...
    'design-%s_', ...
    'model-%s_', ...
    'opts%d.mat'}, '' ...
);

dims = [2, 5, 10];
funs = 1:24;
design = 'ilhs';
inst = 0;

modelResults = struct();
cleanModelResults = struct();

ms = zeros(1, 5);
ma = zeros(1, 5);

for m = 1:nModelTypes
  modelType = S.modelParamDef(m).name;
  nOpts = length(S.modelParamDef(m).values);

  for dim = dims
    for fun = funs
      fun_path = fullfile(res_path, sprintf('%dD', dim), sprintf('f%d', fun));

      for optind = 1:nOpts
        N = 50 * dim;
        fname = sprintf(res_fname_template, dim, fun, inst, N, design, modelType, optind);
        
        fname = fullfile(fun_path, fname);
        if isfile(fname)
          R = load(fname);

          for r = 1:5
            ms(r) = R.results(r).mse;
            ma(r) = R.results(r).mae;
          end
        else
          ms = nan(1, 5);
          ma = nan(1, 5);
        end

        v = S.modelParamDef(m).values{optind};

        if ismember(modelType, {'forest', 'xgb'})
          row = cell2table({dim, fun, optind, N, ms, ma, v.tree_splitFunc, v.tree_splitGainFunc}, ...
                 'VariableNames', {'Dimension','Function','OptionIndex','SampleSize','MSE','MAE','SplitFunction','GainFunction'});

          if isfield(modelResults, modelType)
            modelResults.(modelType) = [modelResults.(modelType); row];
          else
            modelResults.(modelType) = row;
          end
        elseif ismember(modelType, {'gp'})
          row = cell2table({dim, fun, optind, N, ms, ma, v.hypOptions(1).covFcn}, ...
                           'VariableNames', {'Dimension','Function','OptionIndex','SampleSize','MSE','MAE','CovarianceFunction'});

          if isfield(modelResults, modelType)
            modelResults.(modelType) = [modelResults.(modelType); row];
          else
            modelResults.(modelType) = row;
          end
        else
          error('Unknown model type ''%s''', modelType);
        end
      end
    end
  end

  % clean up results
  results = modelResults.(modelType);
  cleanResults = modelResults.(modelType);
  for r=size(results, 1):-1:1
    if all(isnan(results{r, 5}))
      cleanResults(r,:) = [];
    end
  end
  cleanModelResults.(modelType) = cleanResults;
end

Bagging = modelResults.forest;
CleanBagging = cleanModelResults.forest;

Boosting = modelResults.xgb;
CleanBoosting = cleanModelResults.xgb;

Gps = modelResults.gp;
CleanGps = cleanModelResults.gp;

% u = table2array(unique(CleanBagging(:, 3)));
% 
% for m = 1:10
%   BaggingIndex(u(m)) = m;
% end
% 
% u = table2array(unique(XGB(:, 3)));
% 
% for m = 11:15
%   XGBIndex(u(m-10)) = m;
% end
% 
% u = table2array(unique(GP(:, 3)));
% 
% for m = 16:24
%   GPIndex(u(m-15)) = m;
% end