function [paramSignificance, p, tb, stats] = promisingParam(evals, dataSettings, varargin)
% [paramSignificance, p, tb, stats] = promisingParam(evals, dataSettings)
% Identifies promising parameters in settings according to evaluations in
% evals using anovan and multcompare functions.
%
% [...] = promisingParam(evals, dataSettings, settings) identifies 
% parameters using additional function settings. 
%
% Input:
%   evals        - array of evaluation data (output of catEvalSet)
%   dataSettings - array of data settings (output of catEvalSet)
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'Alpha'     - anovan and multcompare significance level
%     'EvalRatio' - fractions of best-data evaluations to be computed
%     'Target'    - target distance to the fitness function optimum
%
% Output:
%   paramSignificance - number of functions, where the parameter value was
%                       significantly best
%   p                 - anovan p-values | cell-array
%   tb                - anovan tables | cell-array
%   stats             - anovan statistics | cell-array
%
% See Also:
%   anovan, catEvalSet, multcompare

  if nargin < 2
    help promisingParam
    return
  end
  
  % parse input
  settings = settings2struct(varargin);
  
  targetValue = defopts(settings, 'Target', 1e-8);
  evalRatio = defopts(settings, 'EvalRatio', 1);
  alpha = defopts(settings, 'Alpha', 0.05);
  
  % important sizes
  [nFunc, nDims, nSettings] = size(evals);
  [nEvals, nInstances] = size(evals{1,1,1});
  nEvalNumbers = length(evalRatio);
  
  % identify settings differences
  [dFields, dValues] = difField(dataSettings);
  nFields = length(dFields);
  % create factors for ANOVA replicating each settings using number of
  % instances
  difFieldVals = cell(1, nFields);
  for fi = 1:nFields
    % numerical factor
    if all(cellfun(@isnumeric, dValues(fi, :)))
      difFieldVals{fi} = repelem([dValues{fi, :}], nInstances);
    % non-numerical factor
    else
      difFieldVals{fi} = cell(1, nSettings*nInstances);
      for s = 1:nSettings
        difFieldVals{fi}((s-1)*nInstances+1 : s*nInstances) = dValues(fi, s);
      end
    end
  end
  
  % initialize ANOVA and comparison cell-arrays
  anovaEvals = cell(1, nEvalNumbers);
  p = cell(nFunc, nDims, nEvalNumbers);
  tb = cell(nFunc, nDims, nEvalNumbers);
  stats = cell(nFunc, nDims, nEvalNumbers);
  paramSignificance = cell(1, nFields);
  for par = 1:nFields
    if isnumeric(dValues{par, 1})
      paramSignificance{par} = zeros(length(unique([dValues{par, :}])), nFunc);
    else
      paramSignificance{par} = zeros(length(unique(dValues(par, :))), nFunc);
    end
  end
  
  % calculate values for comparison in each function, dimension, and 
  % instance
  for f = 1:nFunc
    for d = 1:nDims
      for i = 1:nInstances
        actualEvals = cellfun(@(x) x(:, i), evals(f, d, :), 'UniformOutput', false);
        actualEvals = permute(actualEvals, [1, 3, 2]);
        actualEvals = [actualEvals{:}];
        % number of evaluations needed by the best data to reach the
        % optimum
        minEval = min([nEvals, nEvals - sum(actualEvals < targetValue) + 1]);
        % required numbers of evaluations
        evalNumbers = ceil(evalRatio*minEval);
        % get function values in evalNumbers and save them for each
        % settings and instance
        for e = 1:nEvalNumbers
          anovaEvals{e}(i, :) = actualEvals(evalNumbers(e), :);
        end
      end
      
      % perform anova for each number of evaluations
      for e = 1:nEvalNumbers
        % reshape evaluations to vector
        actualAnovaEvals = reshape(anovaEvals{e}, [1, nInstances*nSettings]);
        % anova input cannot contain NaN value if it does, omit it from
        % significance computation
        if ~any(any(isnan(actualAnovaEvals)))
          % perform anova
          [p{f, d, e}, tb{f, d, e}, stats{f, d, e}] = ...
            anovan(actualAnovaEvals, difFieldVals, ...
                   'Model', 'linear', 'VarNames', dFields, 'Display', 'off');

          % identify significant values
          signifID = p{f, d, e} < alpha;
          if any(signifID)
            signifID = find(signifID');
            for id = signifID
              [mult_c, mult_mat] = multcompare(stats{f, d, e}, ...
                'Dimension', id, 'Display', 'on');
              [~, bestParId] = min(mult_mat(:, 1));
              % add values comparable with the best
              comparableValues = mult_c(mult_c(:, 6) >  alpha ...
                                    & ((mult_c(:, 1) == bestParId) ...
                                    |  (mult_c(:, 2) == bestParId)), 1:2);
              bestParId = [bestParId; comparableValues(comparableValues ~= bestParId)];
              % mark significant values in function f if not all values are
              % significant
              if length(bestParId) < size(mult_mat, 1)
                paramSignificance{id}(bestParId, f) = 1;
                paramSignificance{id}(~ismember(1:size(mult_mat, 1), bestParId), f) = -1;
              end
            end
          end
        end
      end
    end
  end
  
  

end