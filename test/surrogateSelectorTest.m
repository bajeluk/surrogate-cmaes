function tests = surrogateSelectorTest
  tests = functiontests(localfunctions);
end

function testChooseDistantPoints(testCase)

  n = 20;
  zSample = randn(n,2);
  %{
  zSample = [3.5784    1.4090; ...
    2.7694    1.4172; ...
   -1.3499    0.6715; ...
    3.0349   -1.2075; ...
    0.7254    0.7172; ...
   -0.0631    1.6302; ...
    0.7147    0.4889; ...
   -0.2050    1.0347; ...
   -0.1241    0.7269; ...
    1.4897   -0.3034];
  %}
  
  xSample = repmat([1 2], size(zSample,1), 1) + zSample;

  % zSample = repmat([1 1], 10, 1) + randn(4,2);
  xTrain = [ 1.5377 1.3188; ...
    2.8339 -0.3077; ...
   -1.2588  0.5664; ...
    1.8622  1.3426];

  xmean = [1 2]';
  sigma = 1;
  BD = [1 0; 0 1];

  f = figure();
  scatter(xSample(:,1), xSample(:,2), 'r+');
  hold on;
  scatter(xTrain(:,1), xTrain(:,2), 'b+');
  plot(xmean(1), xmean(2), 'ko');

  [xPreSample, zPreSample, xDistant] = SurrogateSelector.chooseDistantPoints(4, xSample, zSample, xTrain, xmean, sigma, BD);

  scatter(xPreSample(:,1), xPreSample(:,2), 'go');
  scatter(xDistant(:,1), xDistant(:,2), 'k+');
end
