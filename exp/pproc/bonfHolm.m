function [corrected_p, h] = bonfHolm(pvalues, alpha)
% bonfHolm Bonferroni-Holm (1979) correction for multiple comparisons.  
% This is a sequentially rejective version of the simple Bonferroni 
% correction for multiple comparisons and strongly controls the family-wise
% error rate at level alpha.
%
% It works as follows:
% 1) All p-values are sorted in order of smallest to largest. m is the
%    number p-values.
% 2) If the 1st p-value is greater than or equal to alpha/m, the procedure
%    is stopped and no p-values are significant.  Otherwise, go on.
% 3) The 1st p-value is declared significant and now the second p-value is
%    compared to alpha/(m-1). If the 2nd p-value is greater than or equal 
%    to alpha/(m-1), the procedure is stopped and no further p-values are 
%    significant.  Otherwise, go on. 
% 4) Et cetera.
%
% As stated by Holm (1979) "Except in trivial non-interesting cases the 
% sequentially rejective Bonferroni test has strictly larger probability of
% rejecting false hypotheses and thus it ought to replace the classical 
% Bonferroni test at all instants where the latter usually is applied."
%
% [corrected_p, h] = bonfHolm(pvalues, alpha)
%
% Input:
%   pvalues - A vector or matrix of p-values. If pvalues is a matrix, it
%             can be of any dimensionality (e.g., 2D, 3D, etc...). | double
%   alpha   - The desired family-wise alpha level (i.e., the probability of
%             rejecting one of more null hypotheses when all null
%             hypotheses are really true). | double scalar | default: 0.05
%
% Output:
%   corrected_p - Bonferroni-Holm adjusted p-values. Any adjusted p-values
%                 less than alpha are significant (i.e., that null
%                 hypothesis is rejected). The adjusted value of the
%                 smallest p-value is p*m. The ith smallest adjusted 
%                 p-value is the max of p(i)*(m-i+1) or adj_p(i-1).
%                 Note, corrected p-values can be greater than 1.
%   h           - A binary vector or matrix of the same dimensionality as
%                 'pvalues'. If the ith element of h is 1, then the ith 
%                 p-value of 'pvalues' is significant.  If the ith element
%                 of 'h' is 0, then the ith p-value of 'pvalues' is NOT 
%                 significant.
%
% Example:
%   >> p = [.56 .22 .001 .04 .01]; %five hypothetical p-values
%   >> [cor_p, h] = bonfHolm(p, .05)
%   cor_p =
%      0.5600    0.4400    0.0050    0.1200    0.0400
%   h =
%       0     0     1     0     1
% 
%   Conclusion: the third and fifth p-values are significant, but not the
%     remaining three.
%
% Reference:
%   Holm, S. (1979) A simple sequentially rejective multiple test
%   procedure. Scandinavian Journal of Statistics. 6, 65-70.
%
% Author:
%   David M. Groppe
%   Kutaslab
%   Dept. of Cognitive Science
%   University of California, San Diego
%   March 24, 2010
%   https://www.mathworks.com/matlabcentral/fileexchange/28303-bonferroni-holm-correction-for-multiple-comparisons

  % check input
  if nargin < 2
    if nargin < 1
      error('You need to provide a vector or matrix of p-values.');
    end
    alpha = 0.05;
  end
  % check input p-values
  if ~isempty(find(pvalues < 0, 1))
    error('Some p-values are less than 0.');
  elseif ~isempty(find(pvalues > 1, 1))
    fprintf('WARNING: Some uncorrected p-values are greater than 1.\n');
  end
  % check input alpha
  assert(alpha > 0 && alpha < 1, 'Alpha must be in (0, 1) interval.')
    
  % sort p-values
  s = size(pvalues);
  if isvector(pvalues)
    if size(pvalues, 1) > 1
     pvalues = pvalues'; 
    end
    [sorted_p, sort_ids] = sort(pvalues);    
  else
    [sorted_p, sort_ids] = sort(reshape(pvalues, 1, prod(s)));
  end
  % indices to return sorted_p to pvalues order
  [~, unsort_ids] = sort(sort_ids);
  % number of tests
  nTests = length(sorted_p);
  mult_fac = nTests:-1:1;
  cor_p_sorted = sorted_p.*mult_fac;
  % Bonferroni-Holm adjusted p-value
  cor_p_sorted(2:nTests) = max([cor_p_sorted(1:nTests-1); cor_p_sorted(2:nTests)]);
  corrected_p = cor_p_sorted(unsort_ids);
  corrected_p = reshape(corrected_p, s);
  % reject null hypothesis when p-values lower than alpha
  h = corrected_p < alpha;
end