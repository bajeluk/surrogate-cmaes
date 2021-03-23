function ft = feature_infocontent(X, y, settings)
% ft = FEATURE_INFOCONTENT(X, y, settings) returns information content of 
% fitness sequences features for dataset [X, y] according to settings.
%
% The Information Content of Fitness Sequences approach quantifies the 
% so-called information content of a continuous landscape, i.e.,
% smoothness, ruggedness, or neutrality. This approach is based on a symbol
% sequence Φ = {φ_1, ..., φ_{n−1} }, with
%           / a, if (y_{i+1} - y_i)/||x_{i+1} - x_i|| < -ε
%     φ_i = - b, if |y_{i+1} - y_i|/||x_{i+1} - x_i|| <= ε
%           \ c, if (y_{i+1} - y_i)/||x_{i+1} - x_i|| >  ε
% This sequence is derived from the objective values y_1 , ..., y_n 
% belonging to the n points x_1 , ..., x_n of a random walk across the 
% landscape and depends on the information sensitivity parameter ε > 0.
% 
% This symbol sequence Φ is aggregated by the information content 
% H(ε) := sum_{i \neq j} p_{ij}*log_6 p_{ij}, where p_{ij} is the 
% probability of having the “block” φ_i φ_j , with φ_i, φ_j \in {a, b, c},
% within the sequence. Note that the base of the logarithm was set to six 
% as this equals the number of possible blocks φ_i φ_j for which φ_i \neq
% φ_j, i.e., φ_i φ_j \in {ab, ba, ac, ca, bc, cb}.
%
% Another aggregation of the information is the so-called partial 
% information content M(ε) := |Φ'|/(n − 1), where Φ' is the symbol sequence
% of alternating a’s and c’s, which is derived from Φ by removing all zeros
% and repeated symbols. These two characteristics are then utilized for 
% computing five features (see Munoz Acosta et al., 2015).
%
% When considering NaN's in y as valid state, the sequence Φ consists of
% φ_i \in {a, b, c, N}, where φ_i = N, if y_{i+1} = NaN or y_i = NaN. Thus,
% H(ε) := sum_{i \neq j} p_{ij}*log_{12} p_{ij} due to increased number of
% possible φ_i φ_j blocks, i.e., φ_i φ_j \in {ab, ba, ac, ca, bc, cb, Na,
% Nb, Nc, aN, bN, cN}.
%
% settings:
%   distance      - distance metric (similar to pdist function) | default:
%                   'euclidean'
%   epsilon       - vector of epsilon values (has to contain 0) | default: 
%                   [0, 10.^linspace(-5, 15, 1000)]
%   nan_state     - consider NaN values in y as valid state | default:
%                   false
%   partial_ratio - parameter r for ratio of partial information
%                   sensitivity | default: 0.5
%   sorting       - sorting strategy for φ calculation | {'nn', 'random'} |
%                   default: 'nn'
%   settling      - parameter s for settling sensitivity | default: 0.05
%
% Features:
%   h_max     - maximum information content: max_ε(H(ε))
%   eps_s     - settling sensitivity: log10(min(ε:H(ε) < s)))
%   eps_max   - argmax_ε(H(ε))
%   m0        - initial partial information: M(ε = 0)
%   eps_ratio - ratio of partial information sensitivity: 
%               log10(max(ε:M(ε) > r * m0))

  if nargin < 3
    if nargin < 2
      help feature_infocontent
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end
  
  % features cannot be computed if the number of points is lower than 2
  if size(X, 1) < 2 || numel(y) < 2
    ft.h_max = NaN;
    ft.eps_s = NaN;
    ft.eps_max = NaN;
    ft.m0 = NaN;
    ft.eps_ratio = NaN;
    return
  end

  % parse settings
  distance = defopts(settings, 'distance', 'euclidean');
  epsilon = defopts(settings, 'epsilon', [0, 10.^linspace(-5, 15, 1000)]);
  epsilon = unique(epsilon);
  if ~ismember(0, epsilon)
    warning('Zero has to be member of epsilon setting. Adding 0 to epsilon.')
    epsilon(end+1) = 0;
  end
  sorting = defopts(settings, 'sorting', 'nn');
  assert(any(strcmp(sorting, {'nn', 'random'})), ...
         '%s sorting setting is not implemented.', sorting)
  settling_param = defopts(settings, 'settling', 0.05);
  eps_ratio_tresh = defopts(settings, 'partial_ratio', 0.5);
  nan_state = defopts(settings, 'nan_state', false);
  
  nData = numel(y);
  nEps = numel(epsilon);
  
  % initialize random number generator to gain identical feature values for
  % identical data
  rng_seed = rng;
  rng(nData, 'twister')
  
  nn_seed = defopts(settings, 'nn_seed', randi(nData));
  
  % TODO: find and process duplicates in data
  % [Z, ix, iz] = unique(X, 'rows');
  % if numel(ix) < nData
  %   X = Z;
  %   for i = 1:numel(ix)
  %     % the following row should be different
  %     z(i) = y(ix(i));
  %   end
  % end
  
  % sort values according to sorting strategy and calculate distances
  if strcmp(sorting, 'nn')
    % nearest neighbors
    [seq, distSeq] = nnSequence(X, nn_seed, distance);
  else
    % random sorting
    seq = randperm(nData)';
    distSeq = arrayfun(@(x) pdist2(X(seq(x), :), X(seq(x+1), :), distance), 1:nData-1)';
  end
  
  % return random number generator settings to original value
  rng(rng_seed)
  
  % for each epsilon calculate information content H and partial
  % information content M
  H = NaN(1, nEps);
  M = NaN(1, nEps);
  for e = 1:nEps
    % compute psi
    psi_val = computePsi(y, seq, distSeq, epsilon(e));
    % compute H, M
    if nan_state
      H(e) = computeNaNH(psi_val);
    else
      H(e) = computeH(psi_val);
    end
    M(e) = computeM(psi_val);
  end
  
  % calculate features
  % maximum information content
  [ft.h_max, h_max_id] = max(H);
  % settling sensitivity
  ft.eps_s = log10(min(epsilon(H < settling_param)));
  % epsilon with maximal H
  ft.eps_max = epsilon(h_max_id);
  % initial partial information
  ft.m0 = M(epsilon == 0);
  % ratio of partial information sensitivity
  ft.eps_ratio = log10(max(epsilon(M > eps_ratio_tresh * ft.m0)));
  
  % ensure features to be non-empty
  ft = repStructVal(ft, @isempty, NaN, 'test');
  
end

function [seq, distSeq] = nnSequence(X, nn_seed, nn_dist)
% calculate nearest-neighbour sequence

  nData = size(X, 1);
  seq = zeros(nData, 1);
  distSeq = zeros(nData - 1, 1);
  
  % calculate distances
  pointDist = squareform(pdist(X, nn_dist));
  pointDist = pointDist + diag(NaN(1, nData));
  
  notUsedId = true(1, nData);
  % initial point
  seq(1) = nn_seed;
  notUsedId(nn_seed) = false;
  notUsedNum = find(notUsedId);
  % find minimal distance point
  [distSeq(1), minId] = min(pointDist(nn_seed, notUsedId));
  % following points
  for p = 2:nData-1
    seq(p) = notUsedNum(minId);
    notUsedId(seq(p)) = false;
    notUsedNum = find(notUsedId);
    % find minimal distance point
    [distSeq(p), minId] = min(pointDist(nn_seed, notUsedId));
  end
  seq(nData) = notUsedNum;
end

function psi_val = computePsi(y, seq, distSeq, epsilon)
% calculate psi values using information sensitivity epsilon
  ratio = diff(y(seq)) ./ distSeq;
  psi_val = NaN(numel(seq) - 1, 1);
  psi_val(ratio < -epsilon) = -1;
  psi_val(abs(ratio) <= epsilon) = 0;
  psi_val(ratio > epsilon) = 1;
end

function H = computeH(psi_val)
% calculate information content ignoring NaN values (original)
  a = psi_val(1:end-1);
  b = psi_val(2:end);
  % calculate probabilities
  p(1) = mean((a == -1) & (b ==  0));
  p(2) = mean((a == -1) & (b ==  1));
  p(3) = mean((a ==  0) & (b == -1));
  p(4) = mean((a ==  0) & (b ==  1));
  p(5) = mean((a ==  1) & (b == -1));
  p(6) = mean((a ==  1) & (b ==  0));
  % non zero probability
  nZp = p~=0;
  % calculate H
  H(nZp) = p(nZp).*log(p(nZp))/log(6);
  H = -sum(H);
end

function H = computeNaNH(psi_val)
% calculate information content taking into account NaN values
  a = psi_val(1:end-1);
  b = psi_val(2:end);
  % calculate probabilities
  p(1)  = mean((a == -1) & (b ==  0));
  p(2)  = mean((a == -1) & (b ==  1));
  p(3)  = mean((a ==  0) & (b == -1));
  p(4)  = mean((a ==  0) & (b ==  1));
  p(5)  = mean((a ==  1) & (b == -1));
  p(6)  = mean((a ==  1) & (b ==  0));
  p(7)  = mean( isnan(a) & (b == -1));
  p(8)  = mean( isnan(a) & (b ==  0));
  p(9)  = mean( isnan(a) & (b ==  1));
  p(10) = mean((a == -1) &  isnan(b));
  p(11) = mean((a ==  0) &  isnan(b));
  p(12) = mean((a ==  1) &  isnan(b));
  % non zero probability
  nZp = p~=0;
  % calculate H
  H(nZp) = p(nZp).*log(p(nZp))/log(numel(p));
  H = -sum(H);
end

function M = computeM(psi_val)
% calculate partial information content
  n = numel(psi_val);
  psi_val = psi_val(psi_val ~= 0);
  psi_val = psi_val([false; diff(psi_val) ~= 0]);
  M = numel(psi_val) / (n - 1);
end