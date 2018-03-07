classdef RandomForestModelTest < ModelTest

  properties (TestParameter)    
    % functional parameters
    fNum = {2}; % {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    dim = {2}; % {2, 5, 10, 20, 40};
    m = {5};
    
    % old model parameters
    
    % modelSpec = {'linear', 'quadratic'};
    % minLeafSize = {10};
    % minGain = {0.1};
    % splitGain0 = {'MSE'};
    % splitGain1 = {'MSE'};
    % growFull = {false, true};
    % fuzziness = {0, 0.1};
    % steepness = {5};
    % boosting = {false, true};
    
    % new model parameters
    
    % forest
    nTrees = {10};
    steepness = {5};
    boosting = {false, true};
    
    % tree
    minGain = {0.2};
    minLeafSize = {10};
    minParentSize = {20};
    maxDepth = {inf};
    growFull = {false, true};
    lossFunc  = {'mse'};
    fuzziness = {0.1};
    
    % predictor
    predictorFunc = {'CombiPolynomial'};
    % predictorFunc = {'Constant', 'LmfitPolynomial', 'Polynomial', ...
    %                  'RegressPolynomial', 'CombiPolynomial'};
    weak_coeff = {NaN}; % ConstantModel
    weak_modelSpec = {{'constant', 'linear', 'purequadratic', 'interactions', 'quadratic'}}; % LmfitPolynomial, Polynomial, RegressPolynomial
    
    % split
    splitFunc = {'Axis'};
    % splitFunc = {'Axis', 'Gaussian', 'HillClimbingOblique', 'KMeans', ...
    %              'PairOblique', 'RandomPolynomial', 'RandomRbf', ...
    %              'ResidualOblique'};
    split_transformationOptions = {struct};
    split_soft = {false};
    split_lambda = {1};
    split_nRepeats = {1}; % RandomSplit
    split_nQuantize = {5}; % AxisSplit, HillClimbingObliqueSplit, PairObliqueSplit
    split_pairFcn = {@(x) x*log(x)}; % PairObliqueSplit
    split_pairRatio = {0.01}; % PairObliqueSplit
    split_discrimType = {{'linear', 'quadratic'}}; % GaussianSplit, KMeansSplit
    split_includeInput = {true}; % GaussianSplit, KMeansSplit
    split_nRandomPerturbations = {10}; % HillClimbingObliqueSplit
    split_kmeans_metric = {'sqeuclidean'}; % KMeansSplit
    split_randrbf_metric = {'euclidean'}; % RandomRbfSplit
    split_degree = {'linear'}; % RandomPolynomialSplit, ResidualObliqueSplit
    
    % splitGain
    splitGain = {'DEMSD'};
    % splitGain = {'DEMSD', 'DENN', 'DE', 'MSE', 'Var'}; % {'DEMSD', 'DENN', 'DE', 'Gradient', 'MSE', 'Var'}; 
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
    function test0(testCase, fNum, dim, m, ...
        minGain, minLeafSize, minParentSize, maxDepth, growFull, lossFunc, fuzziness, ...
        predictorFunc, weak_coeff, weak_modelSpec, ...
        splitFunc, split_transformationOptions, split_soft, split_lambda, ...
        split_nRepeats, split_nQuantize, split_discrimType, split_includeInput, ...
        split_nRandomPerturbations, split_kmeans_metric, split_randrbf_metric, split_degree, ...
        splitGain, splitGain_minSize, splitGain_degree, splitGain_polyMethod, ...
        splitGain_modelFunc, splitGain_weightedGain, splitGain_k, splitGain_regularization, ...
        steepness, boosting, nTrees)
      
      % prepare parameter structure
      params = struct;
      params.modelSpec = '';
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.growFull = growFull;
      params.fuzziness = fuzziness;
      params.steepness = steepness;
      params.boosting = boosting;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      % tree model settings
      forestModelOptions = struct;
      
      % tree options
      forestModelOptions.tree_minGain = minGain;
      forestModelOptions.tree_minLeafSize = minLeafSize;
      forestModelOptions.tree_minParentSize = minParentSize;
      forestModelOptions.tree_maxDepth = maxDepth;
      forestModelOptions.tree_growFull = growFull;
      forestModelOptions.tree_lossFunc = str2func(sprintf('%sLossFunc', lossFunc));
      forestModelOptions.tree_fuzziness = fuzziness;
      
      % weak model options
      forestModelOptions.tree_predictorFunc = str2func(sprintf('%sModel', predictorFunc));
      forestModelOptions.weak_coeff = weak_coeff;
      forestModelOptions.weak_modelSpec = weak_modelSpec;
      
      % split options
      if iscell(splitFunc)
        forestModelOptions.tree_splitFunc = cellfun(@(x) str2func(sprintf('%sSplit', x)), splitFunc, 'UniformOutput', false);
      else
        forestModelOptions.tree_splitFunc = str2func(sprintf('%sSplit', splitFunc));
      end
      forestModelOptions.split_transformationOptions = split_transformationOptions;
      forestModelOptions.split_soft = split_soft;
      forestModelOptions.split_lambda = split_lambda;
      forestModelOptions.split_nRepeats = split_nRepeats;
      forestModelOptions.split_nQuantize = split_nQuantize;
      forestModelOptions.split_discrimType = split_discrimType;
      forestModelOptions.split_includeInput = split_includeInput;
      forestModelOptions.split_nRandomPerturbations = split_nRandomPerturbations;
      forestModelOptions.split_kmeans_metric = split_kmeans_metric;
      forestModelOptions.split_randrbf_metric = split_randrbf_metric;
      forestModelOptions.split_degree = split_degree;
      
      % splitGain options
      forestModelOptions.tree_splitGainFunc = str2func(sprintf('%sSplitGain', splitGain));
      forestModelOptions.splitGain_minSize = splitGain_minSize;
      forestModelOptions.splitGain_degree = splitGain_degree;
      forestModelOptions.splitGain_polyMethod = splitGain_polyMethod;
      forestModelOptions.splitGain_modelFunc = splitGain_modelFunc;
      forestModelOptions.splitGain_weightedGain = splitGain_weightedGain;
      forestModelOptions.splitGain_k = splitGain_k;
      forestModelOptions.splitGain_regularization = splitGain_regularization;      
      
      % forest options
      forestModelOptions.rf_nTrees = nTrees;
      forestModelOptions.rf_boosting = boosting;
      modelFunc = @() RandomForestModel(forestModelOptions);
      
      fprintf('***************** f%02d  %dD  [-%d, %d] *****************\n', fNum, dim, m, m)
      printStructure(params);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, dim, m);
    end
    
  end
end