classdef TreeModel < WeakModel
  properties (Constant, Access = protected)
    nodeTemplate = struct(... % template for nodes
        'leaf', true, ...
        'parent', 0, ...
        'left', 0, ...
        'right', 0, ...
        'depth', 0, ...
        'splitter', [], ...
        'predictor', [] ...
        );
        %'X', [], ...
        %'y', []);
  end
  
  properties %(Access = protected)
    tree_nodes         % list of nodes
    tree_nNodes        % nodes count
    tree_minGain       % minimum gain to split
    tree_minLeafSize   % minimum number of data examples in leaf
    tree_minParentSize % minimum number of data examples in parent
    tree_maxDepth      % maximum depth of the tree
    tree_splitFunc     % cell-array of split functions
    tree_splits        % generators for split functions
    tree_splitGainFunc % function defining splitGain
    tree_splitGain     % evaluator for split functions
    tree_predictorFunc % function which creates a model in leaf
    tree_predictorOpts % options of weak model in leaf
    tree_delPred       % delete predictors in intermediate leaves
    tree_growFull      % grows a full tree then prunes, otherwise prunes during splitting
    tree_lossFunc      % loss function used for pruning
    tree_kfoldPrune    % number of folds used for cross-validated pruning
    tree_fuzziness     % use fuzzy splits (range [0,1])
    tree_inputOptions  % input tree model options
  end
  
  methods
    function obj = TreeModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      obj.tree_minGain = defopts(modelOptions, 'tree_minGain', 1e-2);
      obj.tree_minLeafSize = defopts(modelOptions, 'tree_minLeafSize', 2);
      obj.tree_minParentSize = defopts(modelOptions, 'tree_minParentSize', 10);
      obj.tree_maxDepth = defopts(modelOptions, 'tree_maxDepth', inf);
      obj.tree_growFull = defopts(modelOptions, 'tree_growFull', false);
      obj.tree_lossFunc = defopts(modelOptions, 'tree_lossFunc', @mseLossFunc);
      obj.tree_delPred = defopts(modelOptions, 'tree_delPred', true);
      obj.tree_predictorFunc = defopts(modelOptions, 'tree_predictorFunc', ...
        @ConstantModel);
      % predictor is a weak model of WeakModel class the ancestor of handle
      % class, i.e. it should be created when it is needed and not copied
      obj.tree_predictorOpts = defopts(modelOptions, 'tree_predictorOpts', ...
        modelOptions);
      obj.tree_splitFunc = defopts(modelOptions, 'tree_splitFunc', ...
        {@AxisSplit}); 
      if ~iscell(obj.tree_splitFunc)
        obj.tree_splitFunc = {obj.tree_splitFunc};
      end
      obj.tree_splits = cellfun(@(x) x(modelOptions), obj.tree_splitFunc, 'UniformOutput', false);
      obj.tree_splitGainFunc = defopts(modelOptions, 'tree_splitGainFunc', ...
        @MSESplitGain);
      obj.tree_splitGain = obj.tree_splitGainFunc(modelOptions);
      obj.tree_fuzziness = defopts(modelOptions, 'tree_fuzziness', 0);
      obj.tree_kfoldPrune = defopts(modelOptions, 'tree_kfoldPrune', 5);
      obj.tree_inputOptions = modelOptions;
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)      
      obj.tree_nNodes = 0;
      initialSize = min(2^obj.tree_maxDepth, 2 * round(size(X,1) / obj.tree_minLeafSize));
      obj.tree_nodes = repmat(TreeModel.nodeTemplate, initialSize, 1);
      iNodeRoot = obj.addNode();
      obj.trainModelRecursive(X, y, iNodeRoot, 0);
      % prune fully grown tree
      if obj.tree_growFull
        obj.cvPrune(X, y);
      end
      % delete predictors in internal nodes
      if obj.tree_delPred
        interNodeId = find(~[obj.tree_nodes(:).leaf]);
        for n = interNodeId
          obj.tree_nodes(n).predictor = [];
        end
      end
    end

    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      if obj.tree_fuzziness == 0
        [yPred, sd2] = obj.modelPredictRecursive(X, 1);
      else
        [yPred, sd2] = obj.modelPredictFuzzyRecursive(X, 1);
      end
      if nargout >= 3
        ci = varToConfidence(yPred, sd2);
      end
    end
    
    function N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
        treePredictor = obj.tree_predictorFunc(obj.tree_predictorOpts);
        N = treePredictor.getMinTrainPoints(dim) + 1;
    end
    
    function obj = cvPrune(obj, X, y)
    % cross-validated pruning 
    
      try
        prunedTreeNodes = obj.getPruneLevel(X, y);
      catch err
        warning('TreeModel: Pruning failed with the following error:\n%s', getReport(err))
        prunedTreeNodes = 1 : obj.tree_nNodes;
      end
      % if the tree should be pruned
      if numel(prunedTreeNodes) < obj.tree_nNodes
        inTree = ismember(1:numel(obj.tree_nodes), prunedTreeNodes);
        % replace all pruned nodes by empty nodes
        obj.tree_nodes(~inTree) = ...
          repmat(TreeModel.nodeTemplate, sum(~inTree), 1);
        % find and set up leaves in remaining nodes
        for n = prunedTreeNodes
          % if the child does not know it has a parent => the parent is a
          % leaf now
          if ~obj.tree_nodes(n).leaf && ...
              obj.tree_nodes(obj.tree_nodes(n).left).parent == 0
            obj.tree_nodes(n).leaf = true;
            obj.tree_nodes(n).left = 0;
            obj.tree_nodes(n).right = 0;
          end
          % number of nodes has decreased
          obj.tree_nNodes = numel(prunedTreeNodes);
        end
      end
    end
    
    function obj = prune(obj, X, y)
      obj.pruneRecursive(X, y, 1);
    end
  end
  
  methods (Access = protected)
    function trainModelRecursive(obj, X, y, iNode, depth)
      [n, dim] = size(X);
      % count actual minimal leaf size according to data dimension
      minLeafSize = max(obj.tree_minLeafSize, obj.getMinTrainPoints(dim));
      if depth < obj.tree_maxDepth ...
          && n >= obj.tree_minParentSize ...
          && n >= 2 * minLeafSize ...
          && size(unique(y, 'rows'), 1) >= 2
        best = struct('gain', -inf);
        obj.tree_splitGain = obj.tree_splitGain.updateMinSize(minLeafSize);
        splitGain = obj.tree_splitGain.reset(X, y);
        % TODO: splitGain should return best model here
        % obj.tree_nodes(iNode).predictor = obj.trainPredictor(X, y, SplitGain.splitGain_bestmodelID);
        % find split with maximal gain among different types of splitting
        for iSplit = 1:size(obj.tree_splits, 2)
          split = obj.tree_splits{iSplit}.reset(X, y);
          candidate = split.get(splitGain);
          idx = candidate.splitter(X) <= 0.5;
          if sum(idx) >= minLeafSize && sum(~idx) >= minLeafSize
            if candidate.gain > best.gain
              best = candidate;
            end
          end
        end
        % perform split if possible
        if best.gain > -inf && (best.gain >= obj.tree_minGain || obj.tree_growFull)
          idx = best.splitter(X) <= 0.5;
          
          left = struct('idx', idx, 'iNode', obj.addNode());
          left.X = X(left.idx, :);
          left.y = y(left.idx, :);
          obj.trainModelRecursive(left.X, left.y, left.iNode, depth+1);
          
          right = struct('idx', ~idx, 'iNode', obj.addNode());
          right.X = X(right.idx, :);
          right.y = y(right.idx, :);
          obj.trainModelRecursive(right.X, right.y, right.iNode, depth+1);
          
          % if best.gain < obj.tree_minGain && obj.tree_nodes(left.iNode).leaf && obj.tree_nodes(right.iNode).leaf
            % prune
            % obj.tree_nodes(left.iNode) = obj.nodeTemplate;
            % obj.tree_nodes(right.iNode) = obj.nodeTemplate;
          % else
            % keep it
            obj.tree_nodes(iNode).leaf = false;
            obj.tree_nodes(iNode).splitter = best.splitter;
            obj.tree_nodes(iNode).left = left.iNode;
            obj.tree_nodes(iNode).right = right.iNode;
            
            obj.tree_nodes(left.iNode).parent = iNode;
            % if obj.tree_nodes(left.iNode).leaf
              obj.tree_nodes(left.iNode).predictor = obj.trainPredictor(left.X, left.y, best.leftID);
              obj.tree_nodes(left.iNode).depth = depth + 1;
              %obj.tree_nodes(left.iNode).X = left.X;
              %obj.tree_nodes(left.iNode).y = left.y;
            % end
            
            obj.tree_nodes(right.iNode).parent = iNode;
            % if obj.tree_nodes(right.iNode).leaf
              obj.tree_nodes(right.iNode).predictor = obj.trainPredictor(right.X, right.y, best.rightID);
              obj.tree_nodes(right.iNode).depth = depth + 1;
              %obj.tree_nodes(right.iNode).X = right.X;
              %obj.tree_nodes(right.iNode).y = right.y;
            % end
          % end
        end
      end
      if iNode == 1 % && obj.tree_nodes(iNode).leaf
        % predictor in root
        obj.tree_nodes(iNode).predictor = obj.trainPredictor(X, y, 0);
        %obj.tree_nodes(iNode).X = X;
        %obj.tree_nodes(iNode).y = y;
      end
    end
    
    function predictor = trainPredictor(obj, X, y, modelID)
      predictor = obj.tree_predictorFunc(obj.tree_predictorOpts);
      if isprop(predictor, 'weak_models')
        if modelID == 0
          % choose model with minimal loss
          predictor = predictor.trainModel(X, y);
          y_pred = predictor.modelPredict(X);
          for i = 1:numel(y_pred)
            y_mse(i) = obj.tree_lossFunc(y, y_pred{i});
          end
          [~, modelID] = min(y_mse);
          predictor = predictor.setUseModel(modelID);
          return
        end
        predictor = predictor.setUseModel(modelID);
      end
      predictor = predictor.trainModel(X, y);
    end
    
    function [yPred, sd2] = modelPredictRecursive(obj, X, iNode, pruned)
    % predict objective values recursively
      if nargin < 4
        pruned = [1, find([obj.tree_nodes(:).parent])];
      end
      % return values for leaves or nodes which became leaves through
      % pruning
      if obj.tree_nodes(iNode).leaf || ...
         ~ismember(obj.tree_nodes(iNode).left, pruned)
        [yPred, sd2] = obj.tree_nodes(iNode).predictor.modelPredict(X);
      % or let predict all the children
      else
        yPred = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        idx = obj.tree_nodes(iNode).splitter(X) <= 0.5;
        
        left = struct('idx', idx, 'iNode', obj.tree_nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx), sd2(left.idx)] = obj.modelPredictRecursive(...
                                             X(left.idx, :), left.iNode, pruned);
        end
        
        right = struct('idx', ~idx, 'iNode', obj.tree_nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx), sd2(right.idx)] = obj.modelPredictRecursive(...
                                               X(right.idx, :), right.iNode, pruned);
        end
      end
    end
    
    function nodeErr = nodeErrRecursive(obj, X, y, iNode, nodeErr)
      % predict objective values and computer error recursively in each 
      % node
      
      % node prediction
      [yPred, ~] = obj.tree_nodes(iNode).predictor.modelPredict(X);
      nodeErr(iNode) = obj.tree_lossFunc(y, yPred);
      
      % children prediction
      if ~obj.tree_nodes(iNode).leaf
        idx = obj.tree_nodes(iNode).splitter(X) <= 0.5;
        nodeErr = obj.nodeErrRecursive(X(idx, :), y(idx, :), ...
                                       obj.tree_nodes(iNode).left, ...
                                       nodeErr);
        nodeErr = obj.nodeErrRecursive(X(~idx, :), y(~idx, :), ...
                                       obj.tree_nodes(iNode).right, ...
                                       nodeErr);
      end
      
    end
    
    function [yPred, sd2] = modelPredictFuzzyRecursive(obj, X, iNode)
      if obj.tree_nodes(iNode).leaf
        [yPred, sd2] = obj.tree_nodes(iNode).predictor.modelPredict(X);
      else
        yPred = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        p = obj.tree_nodes(iNode).splitter(X);
        pTresholdLeft = 0.5 - 0.5 * obj.tree_fuzziness;
        pTresholdRight = 0.5 + 0.5 * obj.tree_fuzziness;
        
        left = struct('idx', p <= pTresholdLeft, 'iNode', obj.tree_nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx), sd2(left.idx)] = obj.modelPredictFuzzyRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > pTresholdRight, 'iNode', obj.tree_nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx), sd2(right.idx)] = obj.modelPredictFuzzyRecursive(X(right.idx, :), right.iNode);
        end
        
        both = struct('idx', and(p > pTresholdLeft, p <= pTresholdRight));
        if any(both.idx)
          both.X = X(both.idx, :);
          pRight = p(both.idx);
          pLeft = 1 - pRight;
          [yLeft, sd2Left] = obj.modelPredictFuzzyRecursive(both.X, left.iNode);
          [yRight, sd2Right] = obj.modelPredictFuzzyRecursive(both.X, right.iNode);
          [yPred(both.idx)] = pLeft .* yLeft + pRight .* yRight;
          % TODO can we do this?
          [sd2(both.idx)] = pLeft .* sd2Left + pRight .* sd2Right;
          % [sd2(both.idx)] = pLeft.^2 .* sd2Left + pRight.^2 .* sd2Right;
        end
      end
    end
    
    function pruneRecursive(obj, X, y, iNode)
      if obj.tree_nodes(iNode).left > 0 && obj.tree_nodes(iNode).right > 0 ...
        && obj.tree_nodes(obj.tree_nodes(iNode).left).leaf ...
        && obj.tree_nodes(obj.tree_nodes(iNode).right).leaf
        % internal node having two leaves
        % consider making this node a leaf
        yPred = zeros(size(X, 1), 1);
        p = obj.tree_nodes(iNode).splitter(X);
        pTresholdLeft = 0.5 - 0.5 * obj.tree_fuzziness;
        pTresholdRight = 0.5 + 0.5 * obj.tree_fuzziness;
        
        left = struct('idx', p <= pTresholdLeft, 'iNode', obj.tree_nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > pTresholdRight, 'iNode', obj.tree_nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
        
        both = struct('idx', and(p > pTresholdLeft, p <= pTresholdRight));
        if any(both.idx)
          both.X = X(both.idx, :);
          pRight = p(both.idx);
          pLeft = 1 - pRight;
          [yLeft] = obj.modelPredictFuzzyRecursive(both.X, left.iNode);
          [yRight] = obj.modelPredictFuzzyRecursive(both.X, right.iNode);
          [yPred(both.idx)] = pLeft .* yLeft + pRight .* yRight;
        end
        
        objective = obj.tree_lossFunc(y, yPred);
        
        current = struct;
        current.X = [obj.tree_nodes(left.iNode).X; obj.tree_nodes(right.iNode).X];
        current.y = [obj.tree_nodes(left.iNode).y; obj.tree_nodes(right.iNode).y];
        predictorNew = obj.tree_predictorFunc(obj.tree_predictorOpts);
        predictorNew = predictorNew.trainModel(current.X, current.y);
        objectiveNew = obj.tree_lossFunc(y, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          % remove children
          obj.tree_nodes(left.iNode) = TreeModel.nodeTemplate;
          obj.tree_nodes(right.iNode) = TreeModel.nodeTemplate;
          obj.tree_nodes(iNode).left = 0;
          obj.tree_nodes(iNode).right = 0;
          obj.tree_nodes(iNode).predictor = predictorNew;
          obj.tree_nodes(iNode).X = X;
          obj.tree_nodes(iNode).y = y;
        end
      end
      if obj.tree_nodes(iNode).left == 0 && obj.tree_nodes(iNode).right == 0
        % try to replace complex model in leaf with constant model
        objective = obj.tree_lossFunc(y, obj.tree_nodes(iNode).predictor.modelPredict(X));
        
        predictorNew = ConstantModel(struct);
        predictorNew = predictorNew.trainModel(obj.tree_nodes(iNode).X, obj.tree_nodes(iNode).y);
        objectiveNew = obj.tree_lossFunc(y, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          obj.tree_nodes(iNode).predictor = predictorNew;
        end
      end
    end
    
    function iNode = addNode(obj)
      obj.tree_nNodes = obj.tree_nNodes + 1;
      % generate new nodes if there are no empty available
      if size(obj.tree_nodes, 1) < obj.tree_nNodes
        obj.tree_nodes(obj.tree_nNodes : 2 * obj.tree_nNodes, 1) = ...
          repmat(TreeModel.nodeTemplate, obj.tree_nNodes + 1, 1);
      end
      iNode = obj.tree_nNodes;
    end
    
    function [T, alpha] = getCostComplexity(obj, X, y)
    % get cost-complexity numbers and and appropriate order of the tree 
    % nodes
    
      % calculate error-complexity measure for all nodes
      R = Inf*ones(1, obj.tree_nNodes);
      R = obj.nodeErrRecursive(X, y, 1, R);
    
      % predefine specific function handles
      lt = @(t) obj.tree_nodes(t).left;
      rt = @(t) obj.tree_nodes(t).right;
      
      % init
      alpha_0 = 0; % -realmax;
      for t = obj.tree_nNodes:-1:1
        if obj.tree_nodes(t).leaf
          N(t) = 1;
          S(t) = R(t);
          g(t) = inf;
          G(t) = inf;
        else
          N(t) = N(lt(t)) + N(rt(t));
          S(t) = S(lt(t)) + S(rt(t));
          g(t) = (R(t) - S(t)) / (N(t) - 1);
          G(t) = min([g(t), G(lt(t)), G(rt(t))]);
        end
      end
      
      k = 1;
      T_0 = 1:obj.tree_nNodes;
      while 1 > 0
        % 2.
        if G(1) > alpha_0 + eps
          alpha(k) = alpha_0;
          T{k} = [];
          if k == 1
            T_prev = T_0;
          else
            T_prev = T{k-1};
          end
          for tt = T_prev
            if all(g(obj.getAncestors(tt)) > alpha(k))
              T{k} = [T{k}, tt];
            end
          end
          % uncomment following line for debugging
          % fprintf('%d %d %g %g\n', k, N(1), alpha(k), S(1))
          alpha_0 = G(1);
          k = k + 1;
        end
        if N(1) == 1
          break
        end
        % 3.
        t = 1;
        while G(t) < g(t) - eps
          if G(t) == G(lt(t))
            t = lt(t);
          else
            t = rt(t);
          end
        end
        % 4. Make current node t terminal by setting 
        N(t) = 1;
        S(t) = R(t);
        g(t) = inf;
        G(t) = inf;
        % 5. Update ancestor’s information of current node t
        while t > 1
          t = obj.tree_nodes(t).parent;
          N(t) = N(lt(t)) + N(rt(t));
          S(t) = S(lt(t)) + S(rt(t));
          g(t) = (R(t) - S(t)) / (N(t) - 1);
          G(t) = min([g(t), G(lt(t)), G(rt(t))]);
        end
      end
    end
    
    function [alpha, R] = getLevelCostComplexity(obj, X, y)
    % get cost-complexity numbers of individual pruning levels
      treeDepths = [obj.tree_nodes(1:obj.tree_nNodes).depth];
      depthLevels = unique(treeDepths);
      maxDepth = max(depthLevels);
      R = NaN(maxDepth + 1, 1);
      N = NaN(maxDepth + 1, 1);
      for d = depthLevels
        y_pred = obj.modelPredictRecursive(X, 1, d);
        R(d+1) = obj.tree_lossFunc(y, y_pred);
        % number of leaves per depth level
        N(d+1) = sum(treeDepths == d) + sum([obj.tree_nodes(treeDepths < d).leaf]);
      end
      alpha = (R - R(end))./(N - 1);
    end
    
    function prunedNodes = getPruneLevel(obj, X, y)
    % get cross-validated pruning level
      
      % calculate original tree cost-complexities (alpha)
      [prunedNodes, alpha] = obj.getCostComplexity(X, y(:, 1));
      avgalpha = [sqrt(alpha(1:end-1) .* alpha(2:end))'; Inf];
      nAlpha = numel(avgalpha);
    
      nData = size(X, 1);
      foldId = cvInd(nData, obj.tree_kfoldPrune);
      % change settings not to perform pruning inside cross-validation
      foldTreeOptions = obj.tree_inputOptions;
      foldTreeOptions.tree_growFull = false;
      foldTreeOptions.tree_delPred = false;
      foldTreeOptions.tree_minGain = -inf;
      % create fold tree
      T = TreeModel(foldTreeOptions); 
      
      % prepare objective table
      Y_cv = NaN(nData, nAlpha);
      adjfactor = 1 + 100*eps;
      % cv
      for f = 1:obj.tree_kfoldPrune
        trainId = foldId ~= f;
        testId = ~trainId;
        X_train = X(trainId, :);
        y_train = y(trainId, :);
        X_test = X(testId, :);
        y_test = y(testId, :);
        % train fold tree
        T.trainModel(X_train, y_train);
        % gain cost-complexities
        [T_f, alpha_f] = T.getCostComplexity(X_test, y_test);
        % gain pruning levels
        prunelev = zeros(size(avgalpha));
        for j = 1 : nAlpha
          prunelev(j) = sum(alpha_f <= avgalpha(j)*adjfactor);
        end
        % predict values for chosen pruning levels
        for j = 1 : nAlpha
          Y_cv(testId, j) = T.modelPredictRecursive(X_test, 1, T_f{prunelev(j)});
        end
      end
      
      % calculate cross-validation errors
      err_cv = NaN(nAlpha, 1);
      for j = 1 : nAlpha
        err_cv(j) = obj.tree_lossFunc(y(:, 1), Y_cv(:, j));
      end
      minerr = min(err_cv);
      % return list of nodes after pruning
      prunedNodes = prunedNodes{find(err_cv <= adjfactor*minerr, 1, 'last')};
   
    end
    
    function res = nSubtreeLeaves(obj, t)
    % returns number of current subtree leaves
      if obj.tree_nodes(t).leaf
        res = 1;
      else
        res = obj.nSubtreeLeaves(obj.tree_nodes(t).left) + ...
              obj.nSubtreeLeaves(obj.tree_nodes(t).right);
      end
    end
    
    function ancId = getAncestors(obj, t)
    % returns ancestors of node t
      parent = obj.tree_nodes(t).parent;
      if parent == 0
        ancId = [];
      else
        ancId = [obj.getAncestors(parent), parent];
      end
    end
    
  end
  
end