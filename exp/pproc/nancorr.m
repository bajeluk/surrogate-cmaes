function [coef, corrCoef, reals2, nans2] = nancorr(x, varargin)
%NANCORR Linear or rank correlation taking into account missing values
%   (NaNs) when 'rows' property is set to 'pairwise' as follows:
%     NANRHO(i, j) = R(i, j)*RHO(i, j) + N(i, j), 
%   where R(i, j) is a fraction of rows with no missing values in column i
%   and j, RHO(i, j) is a correlation coefficient between columns i and j,
%   and N(i, j) is a fraction of rows with both missing values in column i
%   and j.
%
%   NANRHO = NANCORR(X) returns a P-by-P matrix containing the pairwise
%   linear correlation coefficient between each pair of columns in the
%   N-by-P matrix X.
%
%   NANRHO = NANCORR(X,Y,...) returns a P1-by-P2 matrix containing the
%   pairwise correlation coefficient between each pair of columns in the 
%   N-by-P1 and N-by-P2 matrices X and Y.
%
%   [NANRHO, RHO, R, N] = NANCORR(X,...) also returns RHO, R, and N
%   matrices. RHO is a matrix of original correlation coefficients, R is a
%   matrix of fractions of rows with no missing values in columns i and j
%   and N is a matrix of fractions of rows, where values are missing in
%   both columns i and j. If parameter 'rows' is not set to 'pairwise', R
%   is a matrix of ones and N is a matrix of zeros.
%
%   [...] = CORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%
%        Parameter  Value
%         'type'    'Pearson' (the default) to compute Pearson's linear
%                   correlation coefficient, 'Kendall' to compute Kendall's
%                   tau, 'Spearman' to compute Spearman's rho, or
%                   'Schweizer' to compute Schweizer-Wolff correlation
%                   coefficient.
%         'rows'    'all' (default) to use all rows regardless of missing
%                   values (NaNs), 'complete' to use only rows with no
%                   missing values, or 'pairwise' to compute RHO(i,j) using
%                   rows with no missing values in column i or j.
%         'tail'    The alternative hypothesis against which to compute
%                   p-values for testing the hypothesis of no correlation.
%                   Choices are:
%                      TAIL         Alternative Hypothesis
%                   ---------------------------------------------------
%                     'both'     correlation is not zero (the default)
%                     'right'    correlation is greater than zero
%                     'left'     correlation is less than zero
%
%   The 'pairwise' option for the 'rows' parameter can produce RHO that is
%   not positive definite.  The 'complete' option always produces a
%   positive definite RHO, but when data are missing, the estimates may be
%   based on fewer observations.
%
%   See also CORR, CORRCOEF, PARTIALCORR, CORRSCHWEIZER

  if nargin < 1 || isempty(x)
    error('scmaes:nancorr:TooFewInputs', 'Requires a data matrix X.');
  end

  % parse input
  p = inputParser;
  
  addRequired(p, 'x', @isnumeric)
  addOptional(p, 'y', [], @isnumeric)
  addParameter(p, 'rows', 'all', @ischar)
  addParameter(p, 'tail', 'both', @ischar)
  addParameter(p, 'type', 'Spearman', @ischar)
  
  parse(p, x, varargin{:});
  
  y = p.Results.y;

  [nRows, nColsX] = size(x);
  [nRowsY, nColsY] = size(y);
  % check y size
  if ~isempty(y)
    assert(nRows == nRowsY, 'scmaes:nancorr:unequalNumOfRows', ...
      'Number of X and Y rows are not equal.')
  else
    nColsY = nColsX;
  end

  % choose according to 'type'
  if strcmpi(p.Results.type, 'schweizer')
    corrCoef = corrSchweizer(x, varargin{:});
  else
    corrCoef = corr(x, varargin{:});
  end

  % calculate NaN coefficients
  reals2 = zeros(nColsX, nColsY);
  nans2 = zeros(nColsX, nColsY);
  xnan = isnan(x);
  ynan = isnan(y);
  switch p.Results.rows
    % calculate numbers of NaN values pairwise
    case 'pairwise'
      % one data matrix input
      if isempty(y)
        for i = 1:nColsX-1
          for j = i+1:nColsX
            reals2(i,j) = sum(~xnan(:, i) & ~xnan(:, j));
            nans2(i,j)  = sum( xnan(:, i) &  xnan(:, j));
          end
        end
        reals2 = reals2 + reals2' + nRows*eye(nColsX);
        nans2 = nans2 + nans2';
      % two data matrices input
      else
        for i = 1:nColsX
          for j = 1:nColsY
            reals2(i,j) = sum(~xnan(:, i) & ~ynan(:, j));
            nans2(i,j)  = sum( xnan(:, i) &  ynan(:, j));
          end
        end
      end
      % calculate coefficients as fractions
      reals2 = reals2/nRows;
      nans2 = nans2/nRows;
    otherwise
      reals2 = ones(nColsX, nColsY);
  end

  % calculate similarity measure which takes into account also the numbers 
  % of NaN values as a sum of following members:
  %   (fraction of both x and y values in R)*corrSW(x, y)
  %   (fraction of both x and y values are NaN)*1
  %   (fraction of x and y values, where one is NaN and the other is in R)*0
  coef = reals2.*corrCoef + nans2;
  
  % replace diagonal by exact ones in case of one matrix input
  if isempty(p.Results.y)
    coef(logical(eye(size(coef)))) = 1;
  end

end