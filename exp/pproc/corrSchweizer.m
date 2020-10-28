function coef = corrSchweizer(x, varargin)
%CORRSCHWEIZER Schweizer-Wolff correlation.
%   RHO = CORRSCHWEIZER(X) returns a P-by-P matrix containing the pairwise
%   Schweizer-Wolff correlation coefficient between each pair of columns in
%   the N-by-P matrix X.
%
%   RHO = CORRSCHWEIZER(X,Y,...) returns a P1-by-P2 matrix containing the
%   pairwise correlation coefficient between each pair of columns in the 
%   N-by-P1 and N-by-P2 matrices X and Y.
%
%   RHO = CORRSCHWEIZER(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values.  Valid parameters are the
%   following:
%
%        Parameter  Value
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
%   See also CORR, CORRCOEF, PARTIALCORR.

%   References:
%      [1] Schweizer, B.; Wolff, E. F. On Nonparametric Measures of
%          Dependence for Random Variables. Ann. Statist. 9 (1981), no. 4,
%          879--885.

  if nargin < 1 || isempty(x)
    error('scmaes:corrSchweizer:TooFewInputs', 'Requires a data matrix X.');
  end

  % find 'type' in varargin
  typeId = cellfun(@(x) strcmp('type', x), varargin);
  % remove 'type' from varargin
  if any(typeId)
    % remove 'type'
    varargin(typeId) = [];
    % remove 'type' value
    varargin(typeId) = [];
  end
  % calculate Spearman correlation coefficient
  spearmanCoef = corr(x, varargin{:}, 'type', 'Spearman');
  % calculate Schweizer-Wolff correlation coefficient
  coef = 6/pi*asin(abs(spearmanCoef)/2);
  % replace diagonal by exact ones in case of one matrix input
  if isempty(varargin) || (numel(varargin) > 0 && ~isnumeric(varargin{1}))
    coef(logical(eye(size(coef)))) = 1;
  end

end