classdef TreeModel < ImplModel
  properties (Constant, Access = private)
    nodeTemplate = struct(... % template for nodes
        'parent', 0, ...
        'left', 0, ...
        'right', 0, ...
        'splitter', [], ...
        'predictor', [], ...
        'X', [], ...
        'y', []);
  end
  
  properties %(Access = protected)
    nodes % list of nodes
    nNodes % nodes count
    minGain % minimum gain to split
    minLeafSize % minimum number of data examples in leaf
    minParentSize % minimum number of data examples in parent
    maxDepth % maximum depth of the tree
    splitGenerators % generators for split functions
    splitEvaluator % evaluator for split functions
    predictorFunc % function which creates a model in leaf
    prune % grows a full tree then prunes, otherwise prunes during splitting
  end
  
  methods
    function obj = TreeModel(modelOptions, xMean)
      % constructor
      obj = obj@ImplModel(modelOptions, xMean);
      
      % specific model options
      obj.minGain = defopts(modelOptions, 'minGain', 1e-3);
      obj.minLeafSize = defopts(modelOptions, 'minLeafSize', 2);
      obj.minParentSize = defopts(modelOptions, 'minParentSize', 10);
      obj.maxDepth = defopts(modelOptions, 'maxDepth', inf);
      obj.prune = defopts(modelOptions, 'prune', false);
      obj.predictorFunc = defopts(modelOptions, 'predictorFunc', ...
        @(xMean) ConstantModel(struct, xMean));
      obj.splitGenerators = defopts(modelOptions, 'splitGenerators', ...
        {AxisTreeSplitGenerator(1, 1, false)});
      obj.splitEvaluator = defopts(modelOptions, 'splitEvaluator', ...
        GenericModelTreeSplitEvaluator(@gainMse, obj.predictorFunc));
    end
    
    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      nData = obj.minLeafSize;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the model based on the data (X,y)
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      
      obj.nNodes = 0;
      initialSize = min(2^obj.maxDepth, 2 * round(size(X,1) / obj.minLeafSize));
      obj.nodes = repmat(TreeModel.nodeTemplate, initialSize, 1);
      iNodeRoot = obj.addNode();
      obj.trainModelRecursive(X, y, iNodeRoot, 0);
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      [y, sd2] = obj.modelPredictRecursive(X, 1);
    end
  end
  
  methods (Access = private)
    function trainModelRecursive(obj, X, y, iNode, depth)
      [N, D] = size(X);
      if depth < obj.maxDepth && N >= obj.minParentSize && N >= 2 * obj.minLeafSize && length(unique(y)) >= 2
        best = struct('gain', -inf);
        obj.splitEvaluator.reset(X, y);
        for iSplitGenerator = 1:size(obj.splitGenerators, 2)
          splitGenerator = obj.splitGenerators{iSplitGenerator};
          splitGenerator.reset(X, y);
          while splitGenerator.hasNext()
            candidate = struct('gain', -inf);
            candidate.splitter = splitGenerator.next();
            idx = candidate.splitter(X);
            if sum(idx) >= obj.minLeafSize && sum(~idx) >= obj.minLeafSize
              candidate.gain = obj.splitEvaluator.eval(candidate.splitter);
              if candidate.gain > best.gain
                best = candidate;
              end
            end
          end
        end
        if best.gain > obj.minGain
          obj.nodes(iNode).splitter = best.splitter;
          idx = best.splitter(X);
          
          left = struct('idx', idx, 'iNode', obj.addNode());
          obj.nodes(iNode).left = left.iNode;
          obj.nodes(left.iNode).parent = iNode;
          obj.trainModelRecursive(X(left.idx, :), y(left.idx), left.iNode, depth+1);
          
          right = struct('idx', ~idx, 'iNode', obj.addNode());
          obj.nodes(iNode).right = right.iNode;
          obj.nodes(right.iNode).parent = iNode;
          obj.trainModelRecursive(X(right.idx, :), y(right.idx), right.iNode, depth+1);
          
          return;
        end
      end
      if isempty(obj.nodes(iNode).predictor)
        obj.nodes(iNode).predictor = obj.predictorFunc(mean(X));
        obj.nodes(iNode).predictor.trainModel(X, y, mean(X), 0);
        obj.nodes(iNode).X = X;
        obj.nodes(iNode).y = y;
      end
    end
    
    function [y, sd2] = modelPredictRecursive(obj, X, iNode)
      if isempty(obj.nodes(iNode).splitter)
        [y, sd2] = obj.nodes(iNode).predictor.modelPredict(X);
      else
        y = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        
        left = struct('idx', obj.nodes(iNode).splitter(X), 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [y(left.idx), sd2(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', ~left.idx, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [y(right.idx), sd2(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
      end
    end
    
    function [y, sd2] = modelPredictIterative(obj, X)
      y = zeros(size(X, 1), 1);
      sd2 = zeros(size(X, 1), 1);
      for i = 1:size(X, 1)
        iNode = 1;
        while ~isempty(obj.nodes(iNode).splitter)
          if obj.nodes(iNode).splitter(X(i, :))
            iNode = obj.nodes(iNode).left;
          else
            iNode = obj.nodes(iNode).right;
          end
        end
        [y(i), sd2(i)] = obj.nodes(iNode).predictor.modelPredict(X(i, :));
      end
    end
    
    function iNode = addNode(obj)
      obj.nNodes = obj.nNodes + 1;
      if size(obj.nodes, 1) < obj.nNodes
        obj.children(2 * obj.nNodes, 1) = TreeModel.nodeTemplate;
      end
      iNode = obj.nNodes;
    end
  end
  
end