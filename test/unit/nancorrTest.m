function tests = nancorrTest
% unit test for nancorr function
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase)
  % test empty input
  
  verifyError(testCase, @nancorr, 'scmaes:nancorr:TooFewInputs')
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
  coef = nancorr(X);
  verifySize(testCase, coef, [nCols, nCols])
  % ones on diagonal
  verifyEqual(testCase, diag(coef), ones(nCols, 1))

  % rows settings 'pairwise'
  coef = nancorr(X, 'rows', 'pairwise');
  verifySize(testCase, coef, [nCols, nCols])
  % ones on diagonal
  verifyEqual(testCase, diag(coef), ones(nCols, 1))
  
  % rows settings 'pairwise' with NaN input
  coef = nancorr(Xnan, 'rows', 'pairwise');
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
  coef = nancorr(X, Y);
  verifySize(testCase, coef, [nColsX, nColsY])

  % rows settings 'pairwise'
  coef = nancorr(X, Y, 'rows', 'pairwise');
  verifySize(testCase, coef, [nColsX, nColsY])

  % rows settings 'pairwise' with NaN input
  coef = nancorr(Xnan, Ynan, 'rows', 'pairwise');
  verifySize(testCase, coef, [nColsX, nColsY])
  
  % one input twice should result the same
  coef1 = nancorr(Xnan, 'rows', 'pairwise');
  coef2 = nancorr(Xnan, Xnan, 'rows', 'pairwise');
  verifyLessThan(testCase, abs(coef1-coef2), 10*eps);
end