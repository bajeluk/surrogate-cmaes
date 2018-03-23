classdef TreeModelTest < ModelTest

  properties (TestParameter)
    % functional parameters
    fNum = {2}; % {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    dim = {2, 5, 10, 20, 40};
    m = {5};
    
    % old model parameters
    
    % minGain = {0.2};
    % minLeafSize = {10};
    % modelSpec = {'linear', 'quadratic'};
    % splitGain0 = {'MSE'};
    % splitGain1 = {'DENN'};
    % fuzziness = {0, 0.1};
    % steepness = {5};
    % growFull = {false, true};
    
    % new model parameters
    
    % tree
    minGain = {0.2};
    minLeafSize = {10};
    minParentSize = {20};
    maxDepth = {inf};
    growFull = {false, true};
    lossFunc  = {'mse'};
    fuzziness = {0.1};
    
    % predictor
    predictorFunc = {'Constant', 'LmfitPolynomial', 'Polynomial', ...
                     'RegressPolynomial', 'CombiPolynomial'};
    weak_coeff = {NaN}; % ConstantModel
    weak_modelSpec = {'constant', 'linear', 'purequadratic', 'interactions', 'quadratic'}; % LmfitPolynomial, Polynomial, RegressPolynomial
    
    % split
    splitFunc = {'Axis', 'Gaussian', 'HillClimbingOblique', 'KMeans', ...
                 'PairOblique', 'RandomPolynomial', 'RandomRbf', ...
                 'ResidualOblique'};
    split_transformationOptions = {struct};
    split_soft = {false};
    split_lambda = {1};
    split_maxHyp = {'10*dim'};
    split_nRepeats = {1000}; % RandomSplit
    split_nQuantize = {5}; % AxisSplit, HillClimbingObliqueSplit, PairObliqueSplit
    split_pairFcn = {@(x) x*log(x)}; % PairObliqueSplit
    split_pairRatio = {0.01}; % PairObliqueSplit
    split_discrimType = {{'linear', 'quadratic'}}; % GaussianSplit, KMeansSplit
    split_includeInput = {true}; % GaussianSplit, KMeansSplit
    split_nRandomPerturbations = {10}; % HillClimbingObliqueSplit
    split_axisPhaseHyp = {'ceil(maxHyp/3)'}; % HillClimbingObliqueSplit
    split_kmeans_metric = {'sqeuclidean'}; % KMeansSplit
    split_randrbf_metric = {'euclidean'}; % RandomRbfSplit
    split_degree = {'linear'}; % RandomPolynomialSplit, ResidualObliqueSplit
    
    % splitGain
    splitGain = {'DEMSD', 'DENN', 'DE', 'MSE', 'Var'}; % {'DEMSD', 'DENN', 'DE', 'Gradient', 'MSE', 'Var'}; 
      % GradientSplitGain useful only in case of second derivatives
    splitGain_minSize = {[]};
    splitGain_degree = {[]};
    splitGain_polyMethod = {''};
    splitGain_modelFunc = {@CombiPolynomialModel};
    splitGain_weightedGain = {true};
    splitGain_k = {1}; % DENNSplitGain
    splitGain_regularization = {0}; % GradientSplitGain
    
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
    
    function test2(testCase, fNum, dim, m, ...
        minGain, minLeafSize, minParentSize, maxDepth, growFull, lossFunc, fuzziness, ...
        predictorFunc, weak_coeff, weak_modelSpec, ...
        splitFunc, split_transformationOptions, split_soft, split_lambda, ...
        split_nRepeats, split_nQuantize, split_discrimType, split_includeInput, ...
        split_nRandomPerturbations, split_kmeans_metric, split_randrbf_metric, ...
        split_degree, split_maxHyp, split_axisPhaseHyp, ...
        splitGain, splitGain_minSize, splitGain_degree, splitGain_polyMethod, ...
        splitGain_modelFunc, splitGain_weightedGain, splitGain_k, splitGain_regularization)
      
      % prepare parameter structure
      params = struct;

      % tree parameters
      params.tree_minGain = minGain;
      params.tree_minLeafSize = minLeafSize;
      params.tree_minParentSize = minParentSize;
      params.tree_maxDepth = maxDepth;
      params.tree_growFull = growFull;
      params.tree_lossFunc = lossFunc;
      params.tree_fuzziness = fuzziness;
      % weak model parameters
      params.tree_predictorFunc = predictorFunc;
      params.weak_coeff = weak_coeff;
      params.weak_modelSpec = weak_modelSpec;
      % splitting parameters
      params.tree_splitFunc = splitFunc;
      params.split_transformationOptions = split_transformationOptions;
      params.split_soft = split_soft;
      params.split_lambda = split_lambda;
      params.split_maxHyp = split_maxHyp;
      params.split_nRepeats = split_nRepeats;
      params.split_nQuantize = split_nQuantize;
      params.split_discrimType = split_discrimType;
      params.split_includeInput = split_includeInput;
      params.split_nRandomPerturbations = split_nRandomPerturbations;
      params.split_axisPhaseHyp = split_axisPhaseHyp;
      params.split_kmeans_metric = split_kmeans_metric;
      params.split_randrbf_metric = split_randrbf_metric;
      params.split_degree = split_degree;
      % splitGain parameters
      params.tree_splitGainFunc = splitGain;
      params.splitGain_minSize = splitGain_minSize;
      params.splitGain_degree = splitGain_degree;
      params.splitGain_polyMethod = splitGain_polyMethod;
      params.splitGain_modelFunc = splitGain_modelFunc;
      params.splitGain_weightedGain = splitGain_weightedGain;
      params.splitGain_k = splitGain_k;
      params.splitGain_regularization = splitGain_regularization;
      
      testCase.reset(params, sprintf('_f%02d_%dD', fNum, dim));
      
      % tree model settings
      treeModelOptions = struct;
      
      % tree options
      treeModelOptions.tree_minGain = minGain;
      treeModelOptions.tree_minLeafSize = minLeafSize;
      treeModelOptions.tree_minParentSize = minParentSize;
      treeModelOptions.tree_maxDepth = maxDepth;
      treeModelOptions.tree_growFull = growFull;
      treeModelOptions.tree_lossFunc = str2func(sprintf('%sLossFunc', lossFunc));
      treeModelOptions.tree_fuzziness = fuzziness;
      % weak model options
      treeModelOptions.tree_predictorFunc = str2func(sprintf('%sModel', predictorFunc));
      treeModelOptions.weak_coeff = weak_coeff;
      treeModelOptions.weak_modelSpec = weak_modelSpec;
      % split options
      if iscell(splitFunc)
        treeModelOptions.tree_splitFunc = cellfun(@(x) str2func(sprintf('%sSplit', x)), splitFunc, 'UniformOutput', false);
      else
        treeModelOptions.tree_splitFunc = str2func(sprintf('%sSplit', splitFunc));
      end
      treeModelOptions.split_transformationOptions = split_transformationOptions;
      treeModelOptions.split_soft = split_soft;
      treeModelOptions.split_lambda = split_lambda;
      treeModelOptions.split_maxHyp = split_maxHyp;
      treeModelOptions.split_nRepeats = split_nRepeats;
      treeModelOptions.split_nQuantize = split_nQuantize;
      treeModelOptions.split_discrimType = split_discrimType;
      treeModelOptions.split_includeInput = split_includeInput;
      treeModelOptions.split_nRandomPerturbations = split_nRandomPerturbations;
      treeModelOptions.split_axisPhaseHyp = split_axisPhaseHyp;
      treeModelOptions.split_kmeans_metric = split_kmeans_metric;
      treeModelOptions.split_randrbf_metric = split_randrbf_metric;
      treeModelOptions.split_degree = split_degree;
      % splitGain options
      treeModelOptions.tree_splitGainFunc = str2func(sprintf('%sSplitGain', splitGain));
      treeModelOptions.splitGain_minSize = splitGain_minSize;
      treeModelOptions.splitGain_degree = splitGain_degree;
      treeModelOptions.splitGain_polyMethod = splitGain_polyMethod;
      treeModelOptions.splitGain_modelFunc = splitGain_modelFunc;
      treeModelOptions.splitGain_weightedGain = splitGain_weightedGain;
      treeModelOptions.splitGain_k = splitGain_k;
      treeModelOptions.splitGain_regularization = splitGain_regularization;
      
      fprintf('***************** f%02d  %dD  [-%d, %d] *****************\n', fNum, dim, m, m)
      printStructure(params);
      
      modelFunc = @() TreeModel(treeModelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, dim, m);
    end
  end
end
