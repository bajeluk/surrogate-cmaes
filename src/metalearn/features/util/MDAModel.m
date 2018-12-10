classdef MDAModel < handle
  %MDA -- Mixture Discriminant Analysis
  %   References:
  %     Hastie & Tibshirani, Discriminant Analysis by Gaussian Mixtures, 1996

  properties
    mu     % cluster means
    sigma  % shared cluster covariance
    L      % cholesky L-factor
    J      % number of groups (classes)
    R      % number of clusters per group
    nIter  % number of iterations of EM-algorithm
    jProb  % group prior
    rProb  % mixing probabilities per group
    d      % input space dimensionality
    debug  % print messages
  end

  methods
    function obj = MDAModel(d, J, R, nIter, jProb)
      % constructor
      assert(d > 0);
      obj.d = d;

      assert(J > 1, 'MDAModel: Number of groups must be larger than one.');
      obj.J = J;

      assert(all(R > 0), 'MDAModel: Number of clusters not positive.');
      obj.R = R;

      assert(nIter > 0);
      obj.nIter = nIter;

      assert(all(numel(jProb) == J));
      assert(abs(sum(jProb) - 1) <= 1e-8);
      obj.jProb = jProb;
      
      obj.rProb = zeros(obj.J, max(obj.R));
      obj.mu = zeros(obj.J, max(obj.R), obj.d);
      obj.sigma = zeros(obj.d, obj.d);
      obj.debug = false;
    end

    function fit(obj, X, g)
      n = size(X, 1);
      obj.d = size(X, 2);

      l = zeros(obj.nIter, 1);

      nIterSmallDiff = 0;
      nIterSmallDiffThres = 3;
      % negLDiffThres = -2e0;
      minLDiff = 1e-3;

      for iter = 1:obj.nIter
        % E-step
        
        pc = zeros(n, max(obj.R), obj.J);
        
        if iter == 1
          % init via k-means clustering
          % TODO: other option -- LVQ
          %       init via faster method
          for j = 1:obj.J
            gj = find(g == j);
            Xj = X(gj, :);
            [cj, muj] = kmeans(Xj, obj.R(j));
            obj.mu(j, :, :) = reshape(muj, 1, obj.R(j), obj.d);
            assert(numel(gj) == numel(cj));

            for k = 1:numel(gj)
              pc(gj(k), cj(k), j) = 1;
            end
          end
          assert(all(sum(sum(pc, 2), 3) == 1));
        else
          for j = 1:obj.J
            Y = obj.mu(j, 1:obj.R(j), :);
            assert(all(size(Y) == [1, obj.R(j), obj.d]));
            Y = reshape(Y, obj.R(j), obj.d);
            % squared mahalanobis distance for all (x, y) pairs
            D = (pdist2(X, Y, 'mahalanobis', obj.sigma)).^2;
            assert(all(size(D) == [n, obj.R(j)]));
            
            A = bsxfun(@times, obj.rProb(j, 1:obj.R(j)), exp(-D/2));
            A = A ./ repmat(sum(A, 2), 1, obj.R(j));
            assert(all(size(A) == [n, obj.R(j)]));

            pc(:, 1:obj.R(j), j) = A;
          end
        end

        % M-step
        rProb1 = zeros(obj.J, max(obj.R));
        mu1 = zeros(obj.J, max(obj.R), obj.d);
        sigma1 = zeros(obj.d, obj.d);

        for j = 1:obj.J
          gj = find(g == j);
          pcj = pc(gj, 1:obj.R(j), j);
          assert(all(size(pcj) == [numel(gj), obj.R(j)]));
          
          b = sum(pcj, 1);
          assert(numel(b) == obj.R(j));
          rProb1(j, 1:obj.R(j)) = b / sum(b); % normalized

          Xj = X(gj, :);
          for r = 1:obj.R(j)
            mu1(j, r, :) = sum((Xj .* repmat(pcj(:, r), 1, obj.d)), 1) / b(r);
          end
        end
        
        for j = 1:obj.J
          gj = find(g == j);
          Xj = X(gj, :);

          for r = 1:obj.R(j)
            % assert(all(size(mu1(j, r, :)) == [1 1 obj.d]));
            mujr = reshape(mu1(j, r, :), 1, obj.d);
            diff = Xj - repmat(mujr, numel(gj), 1);
            presigma = bsxfun(@times, permute(diff, [2, 3, 1]), ...
              permute(diff, [3, 2, 1])); % dim x dim x numel(gj)
            % sum over third matrix dimension (gj elements) - see code 
            % below
            sigma1 = sigma1 + sum( ...
              bsxfun(@times, permute(pc(gj, r, j), [3, 2, 1]), presigma), 3);
            
            % previous version:
            % for k = 1:numel(gj)
              % assert(all(size(Xj(k, :)) == size(mujr)));
              % sigma1 = sigma1 + pc(gj(k), r, j).* (diff(k, :)' * diff(k, :));
            % end
          end
        end
        sigma1 = sigma1 / n;

        % update
        obj.rProb = rProb1;
        obj.mu = mu1;
        obj.sigma = sigma1;
        [obj.L, z] = chol(obj.sigma, 'lower');
        if z
          error('MDAModel/fit: updated sigma not positive-definite.');
        end

        % log lik
        A = obj.pxj(X);
        A = A ./ repmat(sum(A, 2), 1, obj.J);
        ll = zeros(n, 1);
        for k = 1:n
          ll(k) = log(A(k, g(k)));
        end
        l(iter) = sum(ll);

        % stopping control
        if iter > 1
          ldiff = l(iter) - l(iter - 1);
%           if ldiff < negLDiffThres
%             error('MDAModel/fit: Negative log lik difference %.4f', ldiff);
%           end

          if abs(ldiff) < minLDiff
            nIterSmallDiff = nIterSmallDiff + 1;
          else
            nIterSmallDiff = 0;
          end
          if nIterSmallDiff >= nIterSmallDiffThres
            if obj.debug
              fprintf('MDAModel/fit: Small log lik difference in %d consecutive iters, stopping.\n', nIterSmallDiff);
            end
            return;
          end
        else
          ldiff = NaN;
        end

        % debug
        if obj.debug
          fprintf('MDAModel/fit: iter: %d, log lik: %.4f, log lik diff: %.4f\n', ...
            iter, l(iter), ldiff);
        end
      end

      if obj.debug
        fprintf('MDAModel/fit: Maximum number of iters reached: %d, stopping.\n', ...
          obj.nIter);
      end
    end

    function G = pxj(obj, X)
      % P(x|j) -- mixing density from clusters in jth group
      n = size(X, 1);
      G = zeros(n, obj.J);
      for j = 1:obj.J
        % computes for all x
        Y = obj.mu(j, 1:obj.R(j), :);
        assert(all(size(Y) == [1, obj.R(j), obj.d]));
        Y = reshape(Y, obj.R(j), obj.d);
        % squared mahalanobis distance for all (x, y) pairs
        D = (pdist2(X, Y, 'mahalanobis', obj.sigma)).^2;

        assert(all(size(D) == [n, obj.R(j)]));
        A = bsxfun(@times, obj.rProb(j, 1:obj.R(j)), exp(-D/2));
        G(:, j) = sqrt((2*pi)^obj.d * det(obj.sigma)) * sum(A, 2);
      end
    end

    function y = predict(obj, X, how)
      % P(j|x)
      n = size(X, 1);

      G = obj.pxj(X);
      P = bsxfun(@times, obj.jProb, G);
      assert(all(size(P) == [n, obj.J]));

      % normalize
      P = P ./ repmat(sum(P, 2), 1, obj.J);

      if nargin < 3 || strcmp(how, 'map')
        % maximum a posteriori
        [~, y] = max(P, [], 2);
      else
        % full posterior
        y = P;
      end
    end
  end

end

