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
  
  % identify settings differences
  [dFields, dValues] = difField(dataSettings);
  
  % calculate values for comparison in each instance
  [nFunc, nDims, nSettings] = size(evals);
  [nEvals, nInstances] = size(evals{1,1,1});
  
  bestEvals = cell(nFunc, nDims, nSettings);
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
        
        % TODO: get function values in evalNumbers and save them for each
        % settings and instance
        % bestEvals{f, d, s} = actualEvals(evalNumbers);
      end
    end
  end

end