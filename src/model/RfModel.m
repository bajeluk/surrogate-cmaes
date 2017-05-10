classdef RfModel < Model
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was built
    trainSigma            % sigma of the generation when the model was built
    trainBD               % BD of the generation when the model was built
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling
    options

    nTrees                % number of regression trees
    nBestPoints           % number of n best training points ordered correctly by prediction of each tree
    minLeaf               % minimum of points in each leaf
    inputFraction         % fraction of points used in training
    nVarToSample          % number of variables to test during training in each decision split
    forest                % ensemble of regression trees
    ordinalRegression     % indicates usage of ordinal regression
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
  end

  methods
    function obj = RfModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'RfModel (constructor): xMean is not a row-vector.');
      obj.options = modelOptions;
      
      % computed values
      obj.useShift  = defopts(obj.options, 'useShift', false);
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      
      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', false);
      
      % forest options
      obj.nTrees = defopts(modelOptions, 'nTrees', 100);
      obj.minLeaf = defopts(modelOptions, 'minLeaf', 2);
      obj.inputFraction = defopts(modelOptions, 'inputFraction', 1);
      obj.nVarToSample = defopts(modelOptions, 'nVarToSample', ceil(obj.dim/3));
      obj.nBestPoints = defopts(modelOptions, 'nBestPoints', 1);
      obj.ordinalRegression = defopts(modelOptions, 'ordinalRegression', false);
      
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: check correctness of the following expression
      nData = obj.minLeaf * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the GP model based on the data (X,y)

      assert(size(xMean,1) == 1, 'RfModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      obj.forest = {};
      [~,yIdx] = sort(y);
      if obj.ordinalRegression
          yTrain = yIdx(yIdx);
      else
          yTrain = y;
      end
      
      % if we want tree elitism
      if obj.nBestPoints > 0
          nBest = min(obj.nBestPoints, length(y));
          sumGoodTrees = 0; allTrees = 0; actualGoodTrees = 0;
          iter = 0; maxTrees = 100*obj.nTrees;
          
          % train new trees until you have not enough good ones
          while ((sumGoodTrees < obj.nTrees) && allTrees < maxTrees)
              iter = iter + 1;
              newForestSize = max([2, min([ceil((obj.nTrees-actualGoodTrees)*allTrees/actualGoodTrees),...
                  obj.nTrees*2^(iter-1), maxTrees-allTrees])]);
              allTrees = allTrees + newForestSize;
              goodTrees = false(1, newForestSize);
              
              % train forest
              trForest = TreeBagger(newForestSize, X, yTrain, 'method', 'regression', ...
                'MinLeaf', obj.minLeaf, ...
                'FBoot', obj.inputFraction, ...
                'NVarToSample', obj.nVarToSample);
              Trees = trForest.Trees;
            
              if (nBest > 0)
              % find trees with elitism
                  parfor treeNum = 1:newForestSize
                      yPredict = predict(Trees{treeNum},X);
                      yPredict = yPredict(yIdx);
                      if ( issorted(yPredict(1:nBest)) && all(yPredict(nBest+1:end)>=yPredict(nBest)) )
                          goodTrees(treeNum) = 1;
                      end
                  end
              else
              % fill the rest with ordinary ones
                  goodTrees = true(1,newForestSize);
              end
              
              % save trees with elitism
              newGoodTrees = sum(goodTrees);
              obj.forest(end+1:end+newGoodTrees) = Trees(goodTrees);
              sumGoodTrees = sumGoodTrees + newGoodTrees;
              actualGoodTrees = actualGoodTrees + newGoodTrees;
              
              if (maxTrees-allTrees == 0 && nBest > 0)
                 fprintf('Cannot create forest with %d best poits. Trying to find %d remaining trees with %d best points.\n', ...
                  nBest,obj.nTrees-sumGoodTrees,nBest-1);
                  if (nBest == 1)
                      maxTrees = (obj.nTrees-sumGoodTrees) + 1;  
                  else
                      maxTrees = 10*(obj.nTrees-sumGoodTrees) + 1;  
                  end
                  allTrees = 1;
                  nBest = nBest - 1;
                  actualGoodTrees = 0;
                  iter = 1;
              end
          end
          
      % simple forest without elitism
      else
          trForest = TreeBagger(obj.nTrees, X, yTrain, 'method', 'regression',...
                'MinLeaf', obj.minLeaf,...
                'FBoot', obj.inputFraction, ...
                'NVarToSample', obj.nVarToSample);
          obj.forest = trForest.Trees;
      end
        
      % count train MSE
      % trainMSE = mean((y - obj.predict(X)).^2);
      % fprintf('  TreeBagger: train MSE = %f\n', trainMSE);
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        trees = obj.forest;
        yPredictions = NaN(size(X,1),obj.nTrees);

        % each tree prediction
        if isa(trees{1}, 'classregtree')
          % for trees trained by classregtree
          parfor treeNum = 1:obj.nTrees
            yPredictions(:,treeNum) = eval(trees{treeNum},XWithShift);
          end
        else
          parfor treeNum = 1:obj.nTrees
            yPredictions(:,treeNum) = predict(trees{treeNum},XWithShift);
          end
        end

        % averaging results
        y = sum(yPredictions,2)./obj.nTrees + obj.shiftY;
        sd2 = sum((yPredictions - repmat(y,1,obj.nTrees)).^2,2)./obj.nTrees;
      else
        y = []; sd2 = [];
        fprintf(2, 'RfModel.predict(): the model is not yet trained!\n');
      end
    end

  end
  
  methods (Access = protected)
      
    function ensemble = trainForest(obj,X,y)
    % trains forest similar to forest made by TreeBagger, but it is three times slower      
      ensemble = cell(1,obj.nTrees);
      minLeaves = obj.minLeaf;
      IF = obj.inputFraction;
      nData = size(X,1);
          
      parfor n = 1:obj.nTrees
        X4Tree = X;
        Xperm = randi(nData,1,nData);
        Xtrain = X4Tree(Xperm(1:round(nData*IF)),:);
        ensemble{n} = fitrtree(Xtrain,y,'MinLeaf',minLeaves,'NumVariablesToSample',ceil(size(X,2)/3),'MinParentSize',2*minLeaves);
      end
    end
  end

end

