function tests = expUtilTest
  tests = functiontests(localfunctions);
end

function testGetParamIndexVector(testCase)
  % for ID=1 it must always generate vector of ones
  verifyEqual(testCase, getParamIndexVector(1, [7, 5, 3]), [1, 1, 1]);
  verifyEqual(testCase, getParamIndexVector(1, [7, 1, 6, 1, 5, 3]), ones(1,6));
  % for ID=prod(nValues), it must always generate the last values
  verifyEqual(testCase, getParamIndexVector(7*5*3, [7, 5, 3]), [7, 5, 3]);
  verifyEqual(testCase, getParamIndexVector(3*7*5*3, [1, 3, 1, 1, 7, 5, 3]), [1,3,1,1,7,5,3]);
  % some other options
  verifyEqual(testCase, getParamIndexVector(23, [3, 5, 2]), [3, 2, 1]);
  verifyEqual(testCase, getParamIndexVector( 7, [2, 2, 2, 2]), [1, 2, 2, 1]);
end

function testGenerateStructOpts(testCase)
  st = struct('one', 1);
  carray = {'f1', {1, 2}, 'f2', { [10, 11], 1:10 }, 'f3', { st } };
  refStruct = struct();
  refStruct(1).f1 = 1; refStruct(1).f2 = [10, 11]; refStruct(1).f3 = st;
  refStruct(2).f1 = 1; refStruct(2).f2 = 1:10;     refStruct(2).f3 = st;
  refStruct(3).f1 = 2; refStruct(3).f2 = [10, 11]; refStruct(3).f3 = st;
  refStruct(4).f1 = 2; refStruct(4).f2 = 1:10;     refStruct(4).f3 = st;

  verifyEqual(testCase, generateStructOpts(carray), refStruct);
end
