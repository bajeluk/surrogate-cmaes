classdef RfModel < Model
  properties    % derived from abstract class "Model"
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
    dataset             % .X and .y
    shiftMean           % vector of the shift in the X-space
    shiftY = 0;         % shift in the f-space
    options

    nTrees
    nBestPoints         % number of n best training points ordered correctly by prediction of each tree
    minLeaf
    inputFraction
    forest
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
      
      % weird learning
      if obj.nBestPoints > 0
          obj.forest = [];
          [~,yIdx] = sort(y,1,'descend');
          badTrees = 0;
          
          % find one tree predicting points in the right order
          while (isempty(obj.forest) && badTrees < obj.nTrees)
              simpleForest = TreeBagger(1,X,y,'method','regression',... 
                'MinLeaf',obj.minLeaf,...
                'FBoot',obj.inputFraction);
              yPredict = predict(simpleForest,X);
              [~,yPredictIdx] = sort(yPredict,1,'descend');
              if (all( yIdx(1:obj.nBestPoints) == yPredictIdx(1:obj.nBestPoints)))
                  obj.forest = simpleForest;
              else
                  badTrees = badTrees + 1;
              end
          end
          
          % if there was not enough elite trees train normal forest
          if badTrees == obj.nTrees
              obj.forest = TreeBagger(obj.nTrees,X,y,'method','regression',... 
            'MinLeaf',obj.minLeaf,...
            'FBoot',obj.inputFraction);
          else
              % append trees predicting points in the right order to the first
              % one
              while obj.forest.NTrees < obj.nTrees
                  simpleForest = TreeBagger(1,X,y,'method','regression',... 
                    'MinLeaf',obj.minLeaf,...
                    'FBoot',obj.inputFraction);
                  yPredict = predict(simpleForest,X);
                  [~,yPredictIdx] = sort(yPredict,1,'descend');
                  if (all( yIdx(1:obj.nBestPoints) == yPredictIdx(1:obj.nBestPoints)))
                      obj.forest = append(obj.forest,simpleForest);
                  end
              end
          end
              
          
          
      else %simple case no elitism
          obj.forest = TreeBagger(obj.nTrees,X,y,'method','regression',... 
            'MinLeaf',obj.minLeaf,...
            'FBoot',obj.inputFraction);
      end
      
      
      % learning
%       obj.forest = TreeBagger(obj.nTrees,X,y,'method','regression',... 
%             'MinLeaf',obj.minLeaf,...
%             'FBoot',obj.inputFraction);
%       
%       correctTrees = [];
%       if obj.nBestPoints > 0
%           [~,yIdx] = sort(y,1,'descend');
%           for i = 1:obj.nTrees
%               yPredict = predict(obj.Forest.Trees{i},X);
%               [~,yPredictIdx] = sort(yPredict,1,'descend');
%               if (all( yIdx(1:obj.nBestPoints) == yPredictIdx(1:obj.nBestPoints)))
%                   correctTrees = [correctTrees i];
%               end
%           end
%           
%           while length(correctTrees) < obj.nTrees
%               
%           end
%       end
        
      % count train MSE
      trainMSE = mean((y - obj.predict(X)).^2);
      fprintf('  TreeBagger: train MSE = %f\n', trainMSE);
    end

    function [y, dev] = predict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        [y,dev] = predict(obj.forest,X);
        y = y + obj.shiftY;
      else
        y = []; dev = [];
        warning('RfModel.predict(): the model is not yet trained!');
      end
    end

  end

end
