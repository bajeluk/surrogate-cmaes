function tests = generateFeaturesTest
  tests = functiontests(localfunctions);
end

function testConstant(testCase)
  X = [2 3 5; 7 11 13];
  XP = generateFeatures(X, 'constant', false);
  verifyEqual(testCase, XP, ...
    []);
  XP = generateFeatures(X, 'constant', true);
  verifyEqual(testCase, XP, ...
    [1; 1]);
  XP = generateFeatures(X, 'constant', true, true);
  verifyEqual(testCase, XP, ...
    [1; 1]);
end

function testLinear(testCase)
  X = [2 3 5; 7 11 13];
  XP = generateFeatures(X, 'linear', false);
  verifyEqual(testCase, XP, ...
    [2 3 5; 7 11 13]);
  XP = generateFeatures(X, 'linear', true);
  verifyEqual(testCase, XP, ...
    [1 2 3 5; 1 7 11 13]);
  XP = generateFeatures(X, 'linear', true, true);
  verifyEqual(testCase, XP, ...
    [2 3 5 1; 7 11 13 1]);
end

function testPureQuadratic(testCase)
  X = [2 3 5; 7 11 13];
  XP = generateFeatures(X, 'purequadratic', false);
  verifyEqual(testCase, XP, ...
    [2 3 5 2^2 3^2 5^2; 7 11 13 7^2 11^2 13^2]);
  XP = generateFeatures(X, 'purequadratic', true);
  verifyEqual(testCase, XP, ...
    [1 2 3 5 2^2 3^2 5^2; 1 7 11 13 7^2 11^2 13^2]);
  XP = generateFeatures(X, 'purequadratic', true, true);
  verifyEqual(testCase, XP, ...
    [2 3 5 2^2 3^2 5^2 1; 7 11 13 7^2 11^2 13^2 1]);
end

function testInteractions(testCase)
  X = [2 3 5; 7 11 13];
  XP = generateFeatures(X, 'interactions', false);
  verifyEqual(testCase, XP, ...
    [2 3 5 2*3 2*5 3*5; 7 11 13 7*11 7*13 11*13]);
  XP = generateFeatures(X, 'interactions', true);
  verifyEqual(testCase, XP, ...
    [1 2 3 5 2*3 2*5 3*5; 1 7 11 13 7*11 7*13 11*13]);
  XP = generateFeatures(X, 'interactions', true, true);
  verifyEqual(testCase, XP, ...
    [2 3 5 2*3 2*5 3*5 1; 7 11 13 7*11 7*13 11*13 1]);
end

function testQuadratic(testCase)
  X = [2 3 5; 7 11 13];
  XP = generateFeatures(X, 'quadratic', false);
  verifyEqual(testCase, XP, ...
    [2 3 5 2*3 2*5 3*5 2*2 3*3 5*5; 7 11 13 7*11 7*13 11*13 7*7 11*11 13*13]);
  XP = generateFeatures(X, 'quadratic', true);
  verifyEqual(testCase, XP, ...
    [1 2 3 5 2*3 2*5 3*5 2*2 3*3 5*5; 1 7 11 13 7*11 7*13 11*13 7*7 11*11 13*13]);
  XP = generateFeatures(X, 'quadratic', true, true);
  verifyEqual(testCase, XP, ...
    [2 3 5 2*3 2*5 3*5 2*2 3*3 5*5 1; 7 11 13 7*11 7*13 11*13 7*7 11*11 13*13 1]);
end

function testOther(testCase)
  X = [2 3 5; 7 11 13];
  XP = generateFeatures(X, 'other', false);
  verifyEqual(testCase, XP, ...
    []);
  XP = generateFeatures(X, 'other', true);
  verifyEqual(testCase, XP, ...
    [1; 1]);
end