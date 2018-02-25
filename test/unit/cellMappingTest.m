function tests = cellMappingTest
% unit test for cell mapping functions
  tests = functiontests(localfunctions);
end

function testCMCell(testCase)
% testing cell mapping objects

  % empty input
  verifyInstanceOf(testCase, CMCell(), 'CMCell')
  % random input without settings
  X = rand(30, 3);
  y = rand(30, 1);
  cm = CMCell(X, y);
  cmFields = {'X', 'y', 'dim', 'lb', 'ub', 'center', 'minx', 'miny', ...
              'maxx', 'maxy'};
  for f = 1:numel(cmFields)
    verifyNotEmpty(testCase, cm.(cmFields{f}))
  end
  % random input with settings
  settings.lb = [-1, 0, min(X(:, 3))];
  settings.ub = [max(X(:, 1)), 1, 0];
  cm = CMCell(X, y, settings);
  for f = 1:numel(cmFields)
    verifyNotEmpty(testCase, cm.(cmFields{f}))
  end
end

function testCreateCMGrid(testCase)
% testing creating cell mapping grid

  % empty input should not generate error
  verifyNotEmpty(testCase, createCMGrid());
end