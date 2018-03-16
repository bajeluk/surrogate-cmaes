classdef PairObliqueSplit < Split
% PairObliqueSplit for each pair of points constructs a normal vector,
% projects all points onto this vector and uses the projected values as
% tresholds for decision to which side the point belongs

  properties %(Access = protected)
    split_nQuantize % quantization of tresholds 
                    %   0 - all tresholds for one pair
                    %   1, 2, 3, ... number of linearly distributed 
                    %   tresholds per pair
    split_pairFcn   % function to compute the number of pairs according to 
                    % the number of points | function handle
  end

  methods
    function obj = PairObliqueSplit(options)
      obj = obj@Split(options);
      obj.split_nQuantize = defopts(options, 'split_nQuantize', 0);
      obj.split_pairFcn   = defopts(options, 'split_pairFcn', @(x) x*log(x) );
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [n, dim] = size(obj.split_X);
      [nPairs, nTreshPerPair] = obj.getNPairsTresh();
      pair = obj.generatePairs(n, nPairs);
      % pair loop
      for p = 1:nPairs
        % normal vector of the hyperplane
        v = obj.split_X(pair(p, 1), :) - obj.split_X(pair(p, 2), :);
        % project X onto v
        values = (obj.split_X * v')';
        % calculate tresholds
        tresholds = obj.calcTresholds(values, dim, nTreshPerPair(p));
        for treshold = tresholds
          candidate = obj.splitCandidate;
          candidate.splitter = obj.createSplitter(@(X) ...
            X * v' - treshold);
          [candidate.gain, candidate.leftID, candidate.rightID] = splitGain.get(candidate.splitter);
          if candidate.gain > best.gain
            best = candidate;
          end
        end
      end
    end
  end
  
  methods (Access = protected)
    function [nPairs, nTreshPerPair] = getNPairsTresh(obj)
    % calculates the number of pairs and the number of tresholds according
    % to the budget of hyperplanes
      [n, dim] = size(obj.split_X);
      % get prescribed number of pairs and tresholds
      nPairs = max(floor(obj.split_pairFcn(n)), 1);
      nTresh = obj.getNQuant(dim);
      % calculate hyperplane budget
      maxHyp = obj.getMaxHyp(n, dim);
      % too many hyperplanes requires adjustment of treshold and pair
      % numbers
      if nPairs*nTresh > maxHyp
        nPairs = min(maxHyp, nPairs);
        nTreshPerPair = floor(maxHyp/nPairs)*ones(1, nPairs);
        nRemain = mod(maxHyp, nPairs);
        % pairs are chosen at random => ordering is arbitrary
        nTreshPerPair = nTreshPerPair + [ones(1, nRemain), zeros(1, nPairs - nRemain)];
      else
        nTreshPerPair = nTresh*ones(1, nPairs);
      end
    end
  end
  
  methods (Access = private, Static)
    function pair = generatePairs(nData, nPairs)
    % generate pairs sequentially
      numAllPairs = nData*(nData-1)/2;
      % TODO: test efficiency of the following condition properly
      if nPairs > floor(numAllPairs/2)
        % all combinations of pairs
        pair = nchoosek(1:nData, 2);
        pairIds = randperm(numAllPairs, nPairs);
        pair = pair(pairIds);
      else
        % sequentially generate pairs
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
      end
    end
  end
  
end