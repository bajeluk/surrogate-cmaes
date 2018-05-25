function tests = forestModelTest
  tests = functiontests(localfunctions);
end

function testConstructor(testCase)
  % test empty input
  try
    m = ForestModel();
  catch err
    verifyError(testCase, @() error(err.identifier, err.message), 'MATLAB:minrhs')
  end
  % test all variables empty
  modelOpts = [];
  xmean = [];
  try
    m = ForestModel(modelOpts, xmean);
  catch err
    verifyError(testCase, @() error(err.identifier, err.message), 'Forest:xnrw')
  end
  % test model options empty
  modelOpts = [];
  xmean = [2, 3];
  m = ForestModel(modelOpts, xmean);
  verifyClass(testCase, m, 'ForestModel');
  
  % test bagging settings
  modelOpts.forestType = 'bagging';
  xmean = [2, 3];
  m = ForestModel(modelOpts, xmean);
  verifyClass(testCase, m, 'ForestModel');
  verifyClass(testCase, m.forestModel, 'RandomForestModel')
  verifyFalse(testCase, m.forestModel.rf_boosting)
  % test boosting settings (lsboost)
  modelOpts.forestType = 'boosting';
  xmean = [2, 3];
  m = ForestModel(modelOpts, xmean);
  verifyClass(testCase, m, 'ForestModel');
  verifyClass(testCase, m.forestModel, 'RandomForestModel')
  verifyTrue(testCase, m.forestModel.rf_boosting)
  % test lsboost settings
  modelOpts.forestType = 'lsboost';
  xmean = [2, 3];
  m = ForestModel(modelOpts, xmean);
  verifyClass(testCase, m, 'ForestModel');
  verifyClass(testCase, m.forestModel, 'RandomForestModel')
  verifyTrue(testCase, m.forestModel.rf_boosting)
  % test xgboost settings
  modelOpts.forestType = 'xgboost';
  xmean = [2, 3];
  m = ForestModel(modelOpts, xmean);
  verifyClass(testCase, m, 'ForestModel');
  verifyClass(testCase, m.forestModel, 'XGBoostModel')
  verifyTrue(testCase, m.forestModel.rf_boosting)
end

function testMemory(testCase)
% testing amount of memory
  dim = 10;
  nData = 50*dim;
  X = rand(nData, dim);
  y = randn(nData, 1);
  
  trainFrac = 4/5;
  nTrain = ceil(nData*trainFrac);
  pointId = randperm(nData);
  X_train = X(pointId(1:nTrain), :);
  y_train = y(pointId(1:nTrain));
  X_test = X(pointId(nTrain+1:end), :);
  y_test = y(pointId(nTrain+1:end));
  
  % global forest settings
  globalModelOpts.rf_nTrees = 100;
  
  % test bagging settings
  modelOpts = globalModelOpts;
  modelOpts.forestType = 'bagging';
  profile clear
  profile -memory on
  m = ForestModel(modelOpts, X_train(1, :));
  m = m.trainModel(X_train, y_train, X_train(1, :), 1);
  y_pred = m.predict(X_test);
  profile viewer
  verifyClass(testCase, m, 'ForestModel');
  verifyClass(testCase, m.forestModel, 'RandomForestModel')
  verifyFalse(testCase, m.forestModel.rf_boosting)
  
  % test xgboost settings
  modelOpts = globalModelOpts;
  modelOpts.forestType = 'xgboost';  
  profile clear
  profile -memory on
  m = ForestModel(modelOpts, X_train(1, :));
  m = m.trainModel(X_train, y_train, X_train(1, :), 1);
  y_pred = m.predict(X_test);
  profile viewer
  verifyClass(testCase, m, 'ForestModel');
  verifyClass(testCase, m.forestModel, 'XGBoostModel')
  verifyTrue(testCase, m.forestModel.rf_boosting)
end