classdef ObliqueHCTreeSplitGenerator < TreeSplitGenerator
  % obj = AxisTreeSplitGenerator(featuresFraction, valuesFraction)
  % for feature = random unique from features
  %   for value = random unique from values(feature)
  %     yield lambda X(:, feature) <= value
  % setting featuresFraction = 1, valuesFraction = 1 generates every
  % possible split point
  
  properties %(Access = protected)
    samplesFraction % fraction of data values to split by
    samples % sampled examples
    generated
    evaluator
    nRepeats
    iRepeats
  end
  
  methods
    function obj = ObliqueHCTreeSplitGenerator(evaluator, samplesFraction, nRepeats)
      if nargin > 0
        obj.samplesFraction = samplesFraction;
        obj.evaluator = evaluator;
        obj.nRepeats = nRepeats;
      end
    end
    
    function reset(obj, X, y)
      % resets the generator with new data
      obj.X = X;
      obj.y = y;
      obj.samples = datasample(X, ...
        ceil(obj.samplesFraction * size(X, 1)), ...
        'Replace', false);
      obj.iRepeats = 1;
    end
    
    function r = hasNext(obj)
      % whether next split function is available
      r = obj.iRepeats <= obj.nRepeats;
    end
    
    function f = next(obj)
      if obj.iRepeats == 1
        axisGenerator = AxisTreeSplitGenerator(1, 1, false);
        axisGenerator.reset(obj.samples, []);
        best = struct('gain', -inf);
        while axisGenerator.hasNext()
          candidate = struct;
          candidate.splitter = axisGenerator.next();
          candidate.gain = obj.evaluator.eval(candidate.splitter);
          if candidate.gain > best.gain
            best = candidate;
          end
        end
        f = functions(best.splitter);
        ws = f.workspace{1};
        H = [zeros(1, size(X, 2)), ws.treshold];
        H(ws.feature) = 1;
      else
        H = rand(1, size(X, 2) + 1);
      end
      obj.nRepeats = obj.nRepeats + 1;
      X1 = [X ones(size(X, 1), 1)]; 
      J = 10;
      while J > 0
        improvement = true;
        while improvement
          improvement = false;
          for d = 1:size(X, 2)
            V = X1 * H;
            U = (H(:, d) * X1(:, d) - V) ./ X1(:, d);
            H1 = H;
            HBest = 0;
            for u = U
              H1(d) = u;
              % if better
              if true
                HBest = H1;
              end
            end
            if true % HBest better than H
              H = HBest;
            else
              if rand() < 0.5
                H = HBest;
              end
            end
          end
        end
        while J > 0
          H1 = H + 0.01 * rand(1, size(X, 2) + 1);
          % if better
          if true
            H = H1;
            break;
          end
          J = J + 1;
        end
      end
    end
  end
end