function tests = archiveTest
  tests = functiontests(localfunctions);
end

function testSaveArchive(testCase)
  X = (0:0.1:2)';
  y = sin(X);

  a = Archive(size(X,2));
  a = a.save(X, y, 1);          % generation 1
  a = a.save(X+2, sin(X+2), 2); % generation 2
  Xtest = [X; X+2];
  ytest = [y; sin(X+2)];
  
  verifyEqual(testCase, a.X, Xtest);
  verifyEqual(testCase, a.y, ytest);
  verifyEqual(testCase, a.gens, [ones(size(X)); 2*ones(size(X))]);
end