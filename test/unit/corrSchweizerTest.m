function tests = corrSchweizerTest
% unit test for Schweizer-Wolff correlation function
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase)
  % test empty input
  
  verifyError(testCase, @corrSchweizer, 'scmaes:corrSchweizer:TooFewInputs')
end

function testOneMatrixInput(testCase)
  % testing one input
    
  % random input
  nRows = 5;
  nCols = 7;
  X = rand(nRows, nCols);
  Xnan = X;
  Xnan(1, 2) = NaN;
  
  % no extra settings
  coef = corrSchweizer(X);
  verifySize(testCase, coef, [nCols, nCols])
  % ones on diagonal
  verifyEqual(testCase, diag(coef), ones(nCols, 1))

  % rows settings 'pairwise'
  coef = corrSchweizer(X, 'rows', 'pairwise');
  verifySize(testCase, coef, [nCols, nCols])
  % ones on diagonal
  verifyEqual(testCase, diag(coef), ones(nCols, 1))
  
  % rows settings 'pairwise' with NaN input
  coef = corrSchweizer(Xnan, 'rows', 'pairwise');
  verifySize(testCase, coef, [nCols, nCols])
  % ones on diagonal
  verifyEqual(testCase, diag(coef), ones(nCols, 1))

end

function testTwoMatrixInput(testCase)
  % testing two inputs
  
  % random inputs
  nRows = 5;
  nColsX = 6;
  nColsY = 7;
  X = rand(nRows, nColsX);
  Y = rand(nRows, nColsY);
  Xnan = X;
  Xnan(1, 2) = NaN;
  Ynan = Y;
  Ynan(2, 3) = NaN;
  
  % no extra settings
  coef = corrSchweizer(X, Y);
  verifySize(testCase, coef, [nColsX, nColsY])

  % rows settings 'pairwise'
  coef = corrSchweizer(X, Y, 'rows', 'pairwise');
  verifySize(testCase, coef, [nColsX, nColsY])

  % rows settings 'pairwise' with NaN input
  coef = corrSchweizer(Xnan, Ynan, 'rows', 'pairwise');
  verifySize(testCase, coef, [nColsX, nColsY])
  
  % one input twice should result the same
  coef1 = corrSchweizer(Xnan, 'rows', 'pairwise');
  coef2 = corrSchweizer(Xnan, Xnan, 'rows', 'pairwise');
  verifyLessThan(testCase, abs(coef1-coef2), 10*eps);
end