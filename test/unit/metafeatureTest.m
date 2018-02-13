function tests = metafeatureTest
% unit test for metafeature functions
  tests = functiontests(localfunctions);
end

function testELADistribution(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_ela_distribution()));
  % test data
  X = rand(10, 3);
  y = randn(10, 1);
  % output fields
  featFields = {'skewness', 'kurtosis', 'number_of_peaks'};
  returnedFields = fieldnames(feature_ela_distribution(X, y));
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testELALevelset(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_ela_levelset()));
  % test data
  X = rand(30, 3);
  y = randn(30, 1);
  qnt = [10, 25, 50];
  % output fields
  featFields = {'mmce_lda_10', 'mmce_lda_25', 'mmce_lda_50', ...
                'mmce_qda_10', 'mmce_qda_25', 'mmce_qda_50', ... 
                'lda_qda_10',  'lda_qda_25',  'lda_qda_50'};
  ft = feature_ela_levelset(X, y);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
  % test mmce values
  for m = 1:numel(returnedFields)
    if ~isempty(strfind(returnedFields{m}, 'mmce'))
      verifyGreaterThanOrEqual(testCase, ft.(returnedFields{m}), 0)
      verifyLessThanOrEqual(testCase, ft.(returnedFields{m}), 1)
    end
  end
  % test combination values
  for q = qnt
    verifyEqual(testCase, ...
                ft.(sprintf('lda_qda_%d', q)), ...
                ft.(sprintf('mmce_lda_%d', q)) / ft.(sprintf('mmce_qda_%d', q)))
  end
end

function testELAMetamodel(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_ela_metamodel()));
  % test data
  X = rand(30, 3);
  y = randn(30, 1);
  % output fields
  featFields = {'lin_simple', 'lin_w_interact', 'quad_simple', 'quad_w_interact'};
  returnedFields = fieldnames(feature_ela_metamodel(X, y));
  printStructure(feature_ela_metamodel(X, y))
 
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testBasic(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_basic()));
  % test data
  X = rand(30, 3);
  y = randn(30, 1);
  % output fields without settings
  featFields = {'dim', 'observations', 'lower_min', 'lower_max', 'upper_min', ...
                'upper_max', 'objective_min', 'objective_max', 'blocks_min', ...
                'blocks_max', 'cells_total', 'cells_filled', 'minimize_fun'};
  returnedFields = fieldnames(feature_basic(X, y));
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end