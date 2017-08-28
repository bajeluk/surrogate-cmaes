function XP = generateFeatures(X, degree, intercept, back)
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
  switch degree
    case 'constant'
      XP = [];
    case 'linear'
      XP = X;
    case 'purequadratic'
      XP = [X X.^2];
    case 'interactions'
      XP = [X nan(n, d*(d-1)/2)];
      dCur = d+1;
      for d1 = 1:d-1
        for d2 = d1+1:d
          XP(:, dCur) = X(:, d1) .* X(:, d2);
          dCur = dCur + 1;
        end
      end
    case 'quadratic'
      XP = [X nan(n, d*(d-1)/2) X.^2];
      dCur = d+1;
      for d1 = 1:d-1
        for d2 = d1+1:d
          XP(:, dCur) = X(:, d1) .* X(:, d2);
          dCur = dCur + 1;
        end
      end
    otherwise
      XP = [];
  end
  if nargin >= 3 && intercept
    I = ones(n, 1);
    if nargin >= 4 && back
      XP = [XP I];
    else
      XP = [I XP];
    end
  end
end