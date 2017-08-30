function [p] = predProbability(y, yPred, sd2)
  % assume normal distribution and compute prediction probability
  p = normpdf(y-yPred, 0, sqrt(sd2));
end