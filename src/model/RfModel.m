classdef RfModel < Model
  properties    % derived from abstract class "Model"
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
    dataset             % .X and .y
    shiftMean           % vector of the shift in the X-space
    shiftY = 0;         % shift in the f-space
    options

    nTrees              % number of regression trees
    nBestPoints         % number of n best training points ordered correctly by prediction of each tree
    minLeaf             % minimum of points in each leaf
    inputFraction       % fraction of points used in training
    forest              % ensemble of regression trees
    ordinalRegression   % indicates usage of ordinal regression
  end

  methods
    function obj = RfModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'RfModel (constructor): xMean is not a row-vector.');
      obj.options = modelOptions;
      obj.dim     = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY  = 0;

      % this is a MOCK/TEST IMPLEMENTATION!
      obj.nTrees = 100;
      obj.minLeaf = 2;
      obj.inputFraction = 1;
      obj.nBestPoints = 1;
      obj.ordinalRegression = false;
      
      if isfield(modelOptions,'nTrees')
          obj.nTrees = modelOptions.nTrees;
      end
      if isfield(modelOptions,'minLeaf')
          obj.minLeaf = modelOptions.minLeaf;
      end
      if isfield(modelOptions,'inputFraction')
          obj.inputFraction = modelOptions.inputFraction;
      end
      if isfield(modelOptions,'nBestPoints')
          obj.nBestPoints = modelOptions.nBestPoints;
      end
       if isfield(modelOptions,'ordinalRegression')
          obj.ordinalRegression = modelOptions.ordinalRegression;
      end
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: *write this* properly according to dimension and
      %       covariance function set in options
      nData = obj.minLeaf * obj.dim;
    end

    function obj = train(obj, X, y, xMean, generation)
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
          sumGoodTrees = 0; allTrees = 0;
          iter = 0; maxTrees = 100*obj.nTrees;
          
          % train new trees until you have not enough good ones
          while ((sumGoodTrees < obj.nTrees) && allTrees < maxTrees)
              iter = iter + 1;
              newForestSize = min([ceil((obj.nTrees-sumGoodTrees)*allTrees/sumGoodTrees),...
                  obj.nTrees*2^(iter-1),maxTrees-allTrees]);
              allTrees = allTrees + newForestSize;
              goodTrees = false(1,newForestSize);
              
              % train forest
              trainForest = TreeBagger(newForestSize,X,yTrain,'method','regression',...
                'MinLeaf',obj.minLeaf,...
                'FBoot',obj.inputFraction);
              % fprintf('Forest with %d trained\n',newForestSize);
              Trees=trainForest.Trees;
            
              % find trees with elitism
              parfor treeNum = 1:newForestSize
                  yPredict = predict(Trees{treeNum},X);
                  yPredict = yPredict(yIdx);
                  if ( issorted(yPredict(1:nBest)) && all(yPredict(nBest+1:end)>=yPredict(nBest)) )
                      goodTrees(treeNum) = 1;
                  end
              end
              
              % save trees with elitism
              newGoodTrees = sum(goodTrees);
              obj.forest(end+1:end+newGoodTrees) = Trees(goodTrees);
              sumGoodTrees = sumGoodTrees + newGoodTrees;
              % fprintf('%d: %d good trees from %d, remaining %d\n',iter,newGoodTrees, newForestSize,obj.nTrees-sumGoodTrees);
          end
          
          % check if we have all trees we wanted, fill the rest with
          % ordinary ones
          if (allTrees == maxTrees && (sumGoodTrees < obj.nTrees))
              fprintf('Cannot create forest with %d best poits. Appending %d remaining ordinary trees.\n', ...
                  obj.nBestPoints,obj.nTrees-sumGoodTrees);
              trainForest = TreeBagger(obj.nTrees-sumGoodTrees,X,yTrain,'method','regression',...
                'MinLeaf',obj.minLeaf,...
                'FBoot',obj.inputFraction);
              obj.forest(end+1:obj.nTrees) = trainForest.Trees;
          else
              obj.forest = obj.forest(1:obj.nTrees);
          end
          
      % simple forest without elitism
      else
          trainForest = TreeBagger(obj.NTrees,X,yTrain,'method','regression',...
                'MinLeaf',obj.minLeaf,...
                'FBoot',obj.inputFraction);
          obj.forest = trainForest.Trees;
      end
        
      % count train MSE
      % trainMSE = mean((y - obj.predict(X)).^2);
      % fprintf('  TreeBagger: train MSE = %f\n', trainMSE);
    end

    function [y, dev] = predict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        trees = obj.forest;
        yPredictions = NaN(size(X,1),obj.nTrees);
        % each tree prediction
        parfor treeNum = 1:obj.nTrees
            yPredictions(:,treeNum) = predict(trees{treeNum},XWithShift);
        end
        % averaging results
        y = sum(yPredictions,2)./obj.nTrees + obj.shiftY;
        dev = sqrt(sum((yPredictions - repmat(y,1,obj.nTrees)).^2)./obj.nTrees);
      else
        y = []; dev = [];
        warning('RfModel.predict(): the model is not yet trained!');
      end
    end

  end

end

