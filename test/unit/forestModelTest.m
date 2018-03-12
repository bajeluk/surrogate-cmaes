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
