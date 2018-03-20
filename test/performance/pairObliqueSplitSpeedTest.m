% PairObliqueSplit creating pairs speed test

posplit = PairObliqueSplit(struct());

% initial settings
numData = 100;
numPair = 100;
numRepeat = 10;
resultDir = fullfile('exp', 'experiments', 'test', 'performance');
resultName = fullfile(resultDir, 'pairObliqueSplitSpeedTest_res.mat');
[~, ~] = mkdir(resultDir);

% initialize variables
nPairs = zeros(numData, numPair);
nDat = 10*(1:numData);
tnGreaterThanTr = NaN(numData, numPair);

if isfile(resultName)
  S = load(resultName);
  datId = ~ismember(nData, nDat);
else
  for d = 1:numel(nDat)
    nData = nDat(d);
    maxPairs = nData*(nData-1)/2;
    nPairs(d, :) = round(linspace(1, maxPairs, numPair));
    for p = 1:numPair
      % existing measurement
      if (p > 1 && nPairs(d, p) == nPairs(d, p-1)) || ...
         (p > 20 && all(tnGreaterThanTr(d, p-20:p-1) == 0))
        tnGreaterThanTr(d, p) = tnGreaterThanTr(d, p-1);
      else
        
        % measure branch using nchoosek
        tic
        % while is used instead of for cycle to be comparable with another
        % while cycle
        t = 1;
        while (t < 1 + numRepeat) && (toc > 0)          
          % all combinations of pairs
          pair = nchoosek(1:nData, 2);
          numAllPairs = nData*(nData-1)/2;
          pairIds = randperm(numAllPairs, nPairs(d, p));
          pair = pair(pairIds);
          t = t + 1;
        end
        tn = toc;
        
        % measure branch using random filling
        t = 1;
        tic
        while (t < 1 + numRepeat) && (toc < tn)
          % sequentially generate pairs at random
          pairsDone = 0;
          pair = [];
          while pairsDone < nPairs
            newPairs = randi(nData, nPairs - pairsDone, 2);
            % find constant pairs
            constPairs = newPairs(:, 1) == newPairs(:, 2);
            % add existing pairs
            newPairs = [pair; newPairs(~constPairs, :)];
            % exclude redundant pairs
            pair = unique(sort(newPairs, 2), 'rows');
            pairsDone = size(pair, 1);
          end
          t = t + 1;
        end
        tr = toc;
        
        % which algorithm was faster?
        tnGreaterThanTr(d, p) = tn > tr;
      end
      fprintf('nData: %d; nPair: %d | %d \n', nData, nPairs(d, p), tnGreaterThanTr(d, p))
    end
    nData = nDat;
    save(resultName, 'nData', 'nPairs', 'tnGreaterThanTr')
  end
  
end

% load created files
if isfile(resultName)
  load(resultName)
else
  error('The file %s could not be loaded', resultName)
end

% analyze results
% smooth first
smoothMat = NaN(size(tnGreaterThanTr));
bestEdge = zeros(1, numData);
for i = 1:numData
  edgeId = find(~tnGreaterThanTr(i, :), 1, 'first');
  if isempty(edgeId)
    edgeId = numPair + 1;
  end
  if edgeId > bestEdge(i)
    bestEdge(i) = edgeId;
  end
  smoothMat(i, :) = [ones(1, bestEdge(i) - 1), zeros(1, numPair - bestEdge(i) + 1)];
end

nPairFrac = bestEdge/numPair;
for i = 1:numData
  nPairVal(i) = nPairs(i, bestEdge(i));
end
% imshow(smoothMat)
% f = fit(nData', nPairFrac', 'exp2');
% plot(f, nData, nPairFrac)

% artificial data for better fit
extraPoints = 0;
divPoints = 15;
% extraData = [nData, 1, 100*(11:10+extraPoints)];
extraData = nData(divPoints:end);
extraPairFrac = [nPairFrac, eps, nPairFrac(end-1)*ones(1, extraPoints)];
% extraPairVal = [nPairVal, 1, nPairVal(end-1)*ones(1, extraPoints)];
extraPairVal = nPairVal(divPoints:end);
% figure(1)
% f = fit(extraData', extraPairFrac', 'exp2');
% plot(f, extraData, extraPairFrac)
figure(2)
% extraData = log(extraData);
% extraPairFrac = log(extraPairFrac);
% f2 = fit(extraData', extraPairFrac', 'poly1');
% plot(f2, extraData, extraPairFrac)
f2 = fit(extraData', extraPairVal', 'poly2');
plot(f2, extraData, extraPairVal)
% lower numbers of points
figure(3)
extraData = [1, nData(1:divPoints)];
extraPairVal = [1, nPairVal(1:divPoints)];
f3 = fit(extraData', extraPairVal', 'poly2');
plot(f3, extraData, extraPairVal)