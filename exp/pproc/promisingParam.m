function promisingParam(evals, dataSettings)
% Identifies promising parameters in settings according to evaluations in
% evals.
%
% Input:
%   evals        - array of evaluation data (output of catEvalSet)
%   dataSettings - array of data settings (output of catEvalSet)
%
% See Also:
%   catEvalSet

  if nargin < 2
    help promisingParam
    return
  end
  
  targetValue = 1e-8;
  evalRatio = 1;
  alpha = 0.05;
  
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
      paramSignificance{par} = zeros(1, length(unique([dValues{par, :}])));
    else
      paramSignificance{par} = zeros(1, length(unique(dValues(par, :))));
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
        % perform anova
        [p{f, d, e}, tb{f, d, e}, stats{f, d, e}] = ...
          anovan(actualAnovaEvals, difFieldVals, ...
                 'Model', 'linear', 'VarNames', dFields);
               
        % identify significant values
        signifID = p{f, d, e} < alpha;
        if any(signifID)
          signifID = find(signifID);
          for id = signifID
            [mult_c, mult_mat] = multcompare(stats{f, d, e}, 'Dimension', id);
            [~, bestParId] = min(mult_mat(:, 1));
            comparableValues = mult_c(mult_c(:, 6) >  alpha ...
                                  & ((mult_c(:, 1) == bestParId) ...
                                  |  (mult_c(:, 2) == bestParId)), 1:2);
            comparableValues = comparableValues(comparableValues ~= bestParId);
            % increase the number of significant functions
            % TODO: the following row throws an error
            paramSignificance{id}(bestParId) = paramSignificance{id}(bestParId) + 1;
            if ~isempty(comparableValues)
              paramSignificance{id}(comparableValues) = paramSignificance{id}(comparableValues) + 1;
            end
          end
        end
      end
    end
  end
  
  

end