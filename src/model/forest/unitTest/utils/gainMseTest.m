function tests = gainMseTest
  tests = functiontests(localfunctions);
end

function testMseNoGain(testCase)
  current = struct('y', [1 1], 'yPred', [0.5 0.5]);
  left = struct('y', [1], 'yPred', [0.5]);
  right = struct('y', [1], 'yPred', [0.5]);
  gain = gainMse(current, left, right);
  verifyLessThan(testCase, abs(gain-0), 1e-3);
end

function testMseGain(testCase)
  current = struct('y', [1 1], 'yPred', [0.5 0.5]);
  left = struct('y', [1], 'yPred', [1]);
  right = struct('y', [1], 'yPred', [1]);
  gain = gainMse(current, left, right);
  verifyLessThan(testCase, abs(gain-0.25), 1e-3);
end
