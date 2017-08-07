function XP = generateFeatures(X, degree, intercept)
  [N, D] = size(X);
  I = []; % intercept
  if intercept
    I = ones(N, 1);
  end
  switch degree
    case 'constant'
      XP = I;
    case 'linear'
      XP = [I X];
    case 'purequadratic'
      XP = [I X X.^2];
    case 'interactions'
      XP = [I X nan(N, D*(D-1)/2)];
      d = size(I, 2)+D+1;
      for d1 = 1:D
        for d2 = d1+1:D
          XP(:, d) = X(:, d1) .* X(:, d2);
          d = d + 1;
        end
      end
    case 'quadratic'
      XP = [I X nan(N, D*(D+1)/2)];
      d = size(I, 2)+D+1;
      for d1 = 1:D
        for d2 = d1:D
          XP(:, d) = X(:, d1) .* X(:, d2);
          d = d + 1;
        end
      end
    otherwise
      XP = [];
  end
end