function paramIndexVector = getParamIndexVector(id, nValues)
% GETPARAMINDEXVECTOR - get the indices which parameters should be used
% from each field (e.g., of the experiment settings)
%
% paramIndexVector = getParamIndexVector(id, nValues)
%
% Input:
%   id      - id of the parameter settings | positive integer
%   nValues - number of values of individual fields | positive integer
%             vector
%
% Output:
%   paramIndexVector - indices of parameter values in particular fields |
%                      positive integer vector
%
% Example:
%
%   Let's have the following settings:
%
%     opts.field1 = {'a', 'b', 'c'};
%     opts.field2 = {'hi', 'hello', 'ciao', 'ahoy', 'tschus'};
%     opts.field3 = { 1, 2};
%
%   To get the settings for id = 23, we have numbers of individual
%   parameter values nValues = [3 5 2], then it follows:
%
%     nValues = [3 5 2] -- |field1| = 3, |field2| = 5, |field3| = 2
%                          (# of values)
%     paramIndexVector = [3 2 1] -- ID means the third value from the
%                                   'field1', the second value from the
%                                   'field2', and the first value from the
%                                   'field3'
%
% See Also:
%   getParamsFromIndex

  if nargin < 1
    help getParamIndexVector
    return
  end

  % bases of the different orders in the multi-base number, e.g.
  % nValues =   [3     3     1     1     2     3     1     1     5     2];
  % orders =  [540   180    60    60    60    30    10    10    10     2];
  orders = ones(size(nValues));
  for i = 1:(length(nValues))
    orders(i) = prod(nValues(i:end));
  end
  % supply final 1 at the end of the orders vector (for not exceed the array
  % in the following while cycle)
  orders(end+1) = 1;

  % we do not handle numeric 0, so we have to start from (id-1)
  x = id - 1;
  i = 1;

  % result -- vector of integer-parameter-values:
  paramIndexVector = zeros(size(nValues));

  % the cycle is similar to conversion to binary on a paper,
  % start with the highest bit
  while (x >= 0 && i <= length(nValues))
    % (1) divide x with the corresponding 'order'
    div = floor(x / orders(i+1));
    % (2) save that number in the final result in the corresponding order
    paramIndexVector(i) = div + 1;
    % (3) use the remainder of the division for the next iteration
    x = mod(x, orders(i+1));
    % (4) switch to the next (lower) order
    i = i + 1;
  end
end
