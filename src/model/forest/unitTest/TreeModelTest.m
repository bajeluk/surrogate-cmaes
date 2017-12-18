classdef TreeModelTest < ModelTest

  properties (TestParameter)
    % functional parameters
    fNum = {2}; %{1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    m = {10};
    
    % old model parameters
    
    % minGain = {0.2};
    % minLeafSize = {10};
    modelSpec = {'linear', 'quadratic'};
    splitGain0 = {'MSE'};
    splitGain1 = {'DENN'};
    % fuzziness = {0, 0.1};
    steepness = {5};
    % growFull = {false, true};
    
    % new model parameters
    
    % tree
    minGain = {0.2};
    minLeafSize = {10};
    minParentSize = {20};
    maxDepth = {inf};
    growFull = {false, true};
    lossFunc  = {'mse'};
    fuzziness = {0, 0.1};
    
    % predictor
    weakFunc = {'Constant', 'LmfitPolynomial', 'Polynomial', ...
                'RegressPolynomial'};
    weak_coeff = {NaN}; % ConstantModel
    weak_modelSpec = {'linear', 'quadratic'}; % LmfitPolynomial, Polynomial, RegressPolynomial
    
    % split
    splitFunc = {'Axis', 'Gaussian', 'HillClimbingOblique', 'KMeans', ...
                 'PairOblique', 'RandomPolynomial', 'RandomRbf', ...
                 'ResidualOblique'};
    split_transformationOptions = {struct};
    split_soft = {false};
    split_lambda = {1};
    split_nRepeats = {1}; % RandomSplit
    split_nQuantize = {0}; % AxisSplit, HillClimbingObliqueSplit, PairObliqueSplit
    split_discrimType = {{'linear', 'quadratic'}}; % GaussianSplit, KMeansSplit
    split_includeInput = {true}; % GaussianSplit, KMeansSplit
    split_nRandomPerturbations = {10}; % HillClimbingObliqueSplit
    split_kmeans_metric = {'sqeuclidean'}; % KMeansSplit
    split_randrbf_metric = {'euclidean'}; % RandomRbfSplit
    split_degree = {'linear'}; % RandomPolynomialSplit, ResidualObliqueSplit
    
    %TODO: splitGain
    splitGain = {'DEMSD', 'DENN', 'DE', 'Gradient', 'MSE', 'Var'};
    splitGain_minSize = {1};
    splitGain_degree = {[]};
    splitGain_polyMethod = {''};
    splitGain_modelFunc = {[]};
    splitGain_weightedGain = {true};
               
               
    
    
    
    steepness = {5};
  end
  
  methods (TestClassSetup)
    function setupClass(testCase)
      testCase.drawEnabled = false;
    end
  end
  
  methods (Test)
%     function test0(testCase, fNum, m, ...
%         minLeafSize, minGain, splitGain0, growFull, fuzziness, steepness)
%       params = struct;
%       params.modelSpec = '';
%       params.tree_minLeafSize = minLeafSize;
%       params.tree_minGain = minGain;
%       params.tree_splitGain = splitGain0;
%       params.tree_growFull = growFull;
%       params.tree_fuzziness = fuzziness;
%       params.tree_steepness = steepness;
%       testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
%       
%       splitGainOptions = struct;
%       splitGainOptions.splitGain_minSize = minLeafSize;
%       splitGainFunc = str2func(sprintf('%sSplitGain', splitGain0));
%       
%       splits = {};
%       splitOptions = struct;
%       splitOptions.split_soft = fuzziness ~= 0;
%       splitOptions.split_lambda = steepness;
%       splits{1} = AxisSplit(splitOptions);
%       
%       treeModelOptions = struct;
%       treeModelOptions.tree_minGain = minGain;
%       treeModelOptions.tree_splitGain = splitGainFunc(splitGainOptions);
%       treeModelOptions.tree_splits = splits;
%       treeModelOptions.tree_growFull = growFull;
%       treeModelOptions.tree_fuzziness = fuzziness;
%       modelFunc = @() TreeModel(treeModelOptions);
%       
%       [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
%     end
%     
%     function test1(testCase, fNum, m, ...
%         modelSpec, minLeafSize, minGain, splitGain1, growFull, fuzziness, steepness)
%       params = struct;
%       params.tree_modelSpec = modelSpec;
%       params.tree_minLeafSize = minLeafSize;
%       params.tree_minGain = minGain;
%       params.tree_splitGain = splitGain1;
%       params.tree_growFull = growFull;
%       params.tree_fuzziness = fuzziness;
%       params.tree_steepness = steepness;
%       testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
%       
%       predictorOptions = struct;
%       predictorOptions.weak_modelSpec = modelSpec;
%       predictorFunc = @() PolynomialModel(predictorOptions);
%       
%       splitGainOptions = struct;
%       splitGainOptions.splitGain_degree = modelSpec;
%       splitGainOptions.splitGain_minSize = minLeafSize;
%       splitGainFunc = str2func(sprintf('%sSplitGain', splitGain1));
%       
%       splits = {};
%       splitOptions = struct;
%       splitOptions.split_soft = fuzziness ~= 0;
%       splitOptions.split_lambda = steepness;
%       splits{1} = AxisSplit(splitOptions);
%       
%       treeModelOptions = struct;
%       treeModelOptions.tree_predictorFunc = predictorFunc;
%       treeModelOptions.tree_minGain = minGain;
%       treeModelOptions.tree_splitGain = splitGainFunc(splitGainOptions);
%       treeModelOptions.tree_splits = splits;
%       treeModelOptions.tree_growFull = growFull;
%       treeModelOptions.tree_fuzziness = fuzziness;
%       modelFunc = @() TreeModel(treeModelOptions);
%       
%       [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
%     end
    
    function test2(testCase, fNum, m, ...
        modelSpec, minLeafSize, minGain, splitGain1, growFull, fuzziness, steepness)
      params = struct;
      params.modelSpec = modelSpec;
      params.tree_minLeafSize = minLeafSize;
      params.tree_minGain = minGain;
      params.tree_splitGainFunc = splitGain1;
      params.tree_growFull = growFull;
      params.tree_fuzziness = fuzziness;
      params.tree_steepness = steepness;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      % tree model settings
      treeModelOptions = struct;
      
      predictorFunc = @PolynomialModel;
      
      treeModelOptions.splitGain_degree = modelSpec;
      treeModelOptions.splitGain_minSize = minLeafSize;
      treeModelOptions.splitGain_k = 3;
      treeModelOptions.tree_splitGainFunc = str2func(sprintf('%sSplitGain', splitGain1));
      
      splits = {};
      treeModelOptions.split_soft = fuzziness ~= 0;
      treeModelOptions.split_lambda = steepness;
      splits{1} = @AxisSplit;
      
      treeModelOptions.tree_predictorFunc = predictorFunc;
      treeModelOptions.tree_minGain = minGain;
      treeModelOptions.tree_splits = splits;
      treeModelOptions.tree_growFull = growFull;
      treeModelOptions.tree_fuzziness = fuzziness;
      modelFunc = @() TreeModel(treeModelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
  end
end