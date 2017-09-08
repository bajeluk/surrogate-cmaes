load('TreeModelTest.mat');
for i = 1:numel(results)
  if isfield(results(i).params, 'modelSpec')
    break
  end
end

results1 = results(1:i-1);
results1 = reshape(results1, 5, numel(results1)/5)';
results2 = results(i:end);
results2 = reshape(results2, 4, numel(results2)/4)';

draw = false;
%%
n = 5;
best = 4;
minTrain = zeros(1, n);
minTest = zeros(1, n);
distTrain = zeros(1, n);
distTest = zeros(1, n);
for i = 1:size(results1)
  res = results1(i,:);
  params = [res.params];
  [minRMSE, m] = min([res.trainRMSE]);
  minTrain(m) = minTrain(m) + 1;
  distTrain = distTrain + ([res.trainRMSE] - minRMSE) / size(results1, 1);
  if m ~= best && draw
    figure;
    scatter(1:n, [res.trainRMSE]);
    set(gca,'XTick', 1:n);
    set(gca,'XTickLabel', {params.splitGain});
    title(sprintf('%s test RMSE', res(best).name));
  end
  [minRMSE, m] = min([res.testRMSE]);
  minTest(m) = minTest(m) + 1;
  distTest = distTest + ([res.testRMSE] - minRMSE) / size(results1, 1);
  if m ~= best && draw
    figure;
    scatter(1:n, [res.testRMSE]);
    set(gca,'XTick', 1:n);
    set(gca,'XTickLabel', {params.splitGain});
    title(sprintf('%s test RMSE', res(best).name));
  end
end

figure;

subplot(2, 2, 1);
scatter(1:n, minTrain);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('How many times train RMSE was the minimum');

subplot(2, 2, 2);
scatter(1:n, minTest);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('How many times test RMSE was the minimum');

subplot(2, 2, 3);
scatter(1:n, distTrain);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('Avg distance to min train RMSE');
min(distTrain)

subplot(2, 2, 4);
scatter(1:n, distTest);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('Avg distance to min test RMSE');
min(distTest)

print('test1.png', '-dpng');

%%
n = 4;
best = 3;
minTrain = zeros(1, n);
minTest = zeros(1, n);
distTrain = zeros(1, n);
distTest = zeros(1, n);
for i = 1:size(results1)
  res = results2(i,:);
  params = [res.params];
  [minRMSE, m] = min([res.trainRMSE]);
  minTrain(m) = minTrain(m) + 1;
  distTrain = distTrain + ([res.trainRMSE] - minRMSE) / size(results2, 1);
  if m ~= best && draw
    figure;
    scatter(1:n, [res.trainRMSE]);
    set(gca,'XTick', 1:n);
    set(gca,'XTickLabel', {params.splitGain});
    title(sprintf('%s test RMSE', res(best).name));
  end
  [minRMSE, m] = min([res.testRMSE]);
  minTest(m) = minTest(m) + 1;
  distTest = distTest + ([res.testRMSE] - minRMSE) / size(results2, 1);
  if m ~= best && draw
    figure;
    scatter(1:n, [res.testRMSE]);
    set(gca,'XTick', 1:n);
    set(gca,'XTickLabel', {params.splitGain});
    title(sprintf('%s test RMSE', res(best).name));
  end
end

figure;

subplot(2, 2, 1);
scatter(1:n, minTrain);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('How many times train RMSE was the minimum');

subplot(2, 2, 2);
scatter(1:n, minTest);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('How many times test RMSE was the minimum');

subplot(2, 2, 3);
scatter(1:n, distTrain);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('Avg distance to min train RMSE');
min(distTrain)

subplot(2, 2, 4);
scatter(1:n, distTest);
set(gca,'XTick', 1:n);
set(gca,'XTickLabel', {params.splitGain});
title('Avg distance to min test RMSE');
min(distTest)

print('test2.png', '-dpng');

%%
for test = 1:2
  minTrain = zeros(1, n);
  minTest = zeros(1, n);
  for fNum = 1:24
    prefix = sprintf('test%d_%02d_010', test, fNum);
    minTrainPos = -1;
    minTestPos = -1;
    for i = 1:numel(results)
      result = results(i);
      if startsWith(result.name, prefix)
        if minTrainPos == -1 || result.trainRMSE < results(minTrainPos).trainRMSE
          minTrainPos = i;
        end
        if minTestPos == -1 || result.testRMSE < results(minTestPos).testRMSE
          minTestPos = i;
        end
      end
    end
    {prefix
      results(minTrainPos).params.splitGain
      results(minTestPos).params.splitGain}
  end
end