function data = divSmooth(data, funcSet, maxFE)
% Divide evaluations by dimension and make data smoother.
%
% Input:
%   data    - cell array functions x dims x nSettings containing cell array
%             of individual instances results
%   funcSet - structure with fields 'BBfunc' (numbers of BBOB
%             functions) and 'dims' (numbers of dimensions) 
%   maxFE   - maximal number of function evaluations
%
% Output:
%   data - cell array of f-values double arrays 
%            - for each function contains maxFE x nI double array with
%              function values, where nI is the number of instances
%
% See Also:
%   dataReady, bbobDataReady

  if nargin < 3
    if nargin < 1
      help divSmooth
      if nargout > 0
        data = [];
      end
      return
    end
    maxFE = 250;
  end

  [func, dims, nSettings] = size(data);
  for s = 1:nSettings
    for d = 1:dims
      for f = 1:func
        if ~isempty(data{f,d,s})
          fInstant = [];
          instances = find(~cellfun(@isempty, data{f,d,s}));
          for i = 1:length(instances)
            data{f,d,s}{instances(i)}(:,2) = ceil(data{f,d,s}{instances(i)}(:,2)/funcSet.dims(d));
            % use only first two columns - fitness, evaluations
            fInstant(:, end+1) = smoothYEvals(data{f,d,s}{instances(i)}(:, 1:2), maxFE);
          end
          data{f,d,s} = fInstant;
        end
      end
    end
  end
  
end