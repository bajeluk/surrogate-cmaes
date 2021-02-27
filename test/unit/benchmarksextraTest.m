function tests = benchmarksextraTest
  tests = functiontests(localfunctions);
end

function testInitAndFinal(testCase)
% test initialization and finalization
  func = 201:207;
  inst = 1;

  for f = func
    opt = fgeneric('initialize', f, inst, '/tmp/deleteme');
    verifySize(testCase, opt, [1, 1])
    fgeneric('finalize');
  end
end

function testF204(testCase)
% test all instances of f204
  inst = 1:10;
  fun = 204;

  % test all instances
  for i = inst
    opt = fgeneric('initialize', fun, i, '/tmp/deleteme');
    verifySize(testCase, opt, [1, 1])
    fgeneric('finalize');
  end
end