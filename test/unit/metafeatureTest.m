function tests = metafeatureTest
% unit test for metafeature functions
  tests = functiontests(localfunctions);
end

function testELADistribution(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_ela_distribution()));
end
