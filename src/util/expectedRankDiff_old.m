function perm = expectedRankDiff_old(yPredict, sd2Predict, mu, varargin)
%EXPECTEDRANKDIFF returns permutation of points which causes highest expected difference in ranking
%
% y       mean predicions returned by the GP/rforest model
% sd2     predicited variances for respective predictions in y
% mu      CMA-ES' mu parameter -- how many points are taken for
%         updating mean and rank-mu update of covariance matrix
%
% perm    Returns permutation of the points in the order from maximum
%         expected error (in perm(1)) to the lowest error (in
%         perm(end)

  if (nargin >= 5)
    rankFunc = varargin{1};
  else
    rankFunc = @errRankMuOnly;
  end

  % sort provided y-values
  [ySort, yInd] = sort(yPredict);
  % sort provided variances, too
  sd2Sort = sd2Predict(yInd);
  % initializations
  n = length(ySort);
  expectedError = zeros(1,n);
  probMatrix = zeros(n,n);

  % calculate expected rank difference for each point
  for i = 1:n
    y = ySort(i);
    sd2 = sd2Sort(i);
    probs = zeros(1,n);

    % probabilities of reaching the values of the other points
    % prob(2:end) = normcdf(ySort([1:(i-1) (i+1):n]), y, sd2);
    probY = 1-normcdf(ySort, y, sd2);
    % probability of being first (unless previously first --
    % something is then added to this value)
    probs(1) = 1-probY(1);
    % probY from the last iteration of the next cycle
    lastProb = probY(1);

    % probabilities of other positions
    position = 1;
    for j = 1:n
      probs(position) = probs(position) + lastProb - probY(j);
      lastProb = probY(j);
      if (i ~= j)
        % find out what is the current inspected permutation
        permutation = 1:n;
        permutation(i) = [];
        permutation = [permutation(1:j-1) i permutation(j:end)];
        % save this contribution to the final expected error
        % of this (i)-th point
        expectedError(i) = expectedError(i) + probs(position) * rankFunc(permutation, mu);

        % we will move on to calculation of the next probability
        position = position + 1;
      end
    end
    % fill also the last probability
    probs(position) = probs(position) + probY(end);
    % save the probabilities into the final probability matrix
    probMatrix(i,:) = probs;

    % Debug:
    % fprintf('Expected error of [%d] is %f.\n', i, expectedError(i));
  end

  % Debug:
  % fprintf('Ranking of errors of sorted points: %s\n', num2str(ranking(-expectedError)));
  % fprintf('Final permutation of original points: %s\n', num2str(perm'));

  % return the final permutation of the points from the highest expected error
  % to the lowest (according to (-1)*expectedError sorted according to
  % inverse sort defined by yInd, which is eqal to ranking(yPredict)
  % yRnk = ranking(yPredict);
  % eRnk = ranking(-expectedError);

  [~, eInd] = sort(-expectedError);
  % the final order of points to reevaluate is following
  perm = yInd(eInd);

  % Debug:
  yRnk = ranking(yPredict);
  assert(all(yRnk(perm)' == eInd), 'Ranking of ''perm'' and ''expectedError'' is not same!');
end

