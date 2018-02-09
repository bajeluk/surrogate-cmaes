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
  X = rand(5, 3);
  y = randn(5, 1);
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