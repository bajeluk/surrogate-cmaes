function XP = generateFeatures(X, degree, intercept)
% Generates new features(columns).
% XP = generateFeatures(X, degree, intercept)
% Generates new features(columns) up to degree and adds intercept.
%
% Input:
%   X           - NxD matrix
%   degree      - polynomial degree
%       'constant'          Intercept only
%       'linear'            Linear (main effect) terms only.
%       'interactions'      Linear and pairwise interaction terms.
%       'purequadratic'     Linear and squared terms.
%       'quadratic'         Linear, pairwise interactions, and squares.
%   intercept   - bool whether to add intercept column (ones)
%
% Output:
%   XP  - NxP matrix where P is the total number of polynomial features
% See Also:
%   fitlm

  [n, d] = size(X);
  I = []; % intercept
  if intercept
    I = ones(n, 1);
  end
  switch degree
    case 'constant'
      XP = I;
    case 'linear'
      XP = [X I];
    case 'purequadratic'
      XP = [X X.^2 I];
    case 'interactions'
      XP = [X nan(n, d*(d-1)/2) I];
      dCur = size(I, 2)+d+1;
      for d1 = 1:d-1
        for d2 = d1+1:d
          XP(:, dCur) = X(:, d1) .* X(:, d2);
          dCur = dCur + 1;
        end
      end
    case 'quadratic'
      XP = [X nan(n, d*(d+1)/2) I];
      dCur = size(I, 2)+d+1;
      for d1 = 1:d
        for d2 = d1:d
          XP(:, dCur) = X(:, d1) .* X(:, d2);
          dCur = dCur + 1;
        end
      end
    otherwise
      XP = [];
  end
end