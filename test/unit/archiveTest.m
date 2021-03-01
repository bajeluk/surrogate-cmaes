function tests = archiveTest
  tests = functiontests(localfunctions);
end

function testSaveArchive(testCase)
% test archive saving ability
  X = (0:0.1:1.9)';
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

function testGetFullData(testCase)
% test getting all data in the archive
  X = randn(20, 2);
  y = fsphere(X);

  a = Archive(size(X, 2));
  a = a.save(X, y, 1);

  % call getFullData function
  [Xtest, yTest] = a.getFullData();

  verifyEqual(testCase, Xtest, X);
  verifyEqual(testCase, yTest, y);
end

function testGetDataNearPoint(testCase)
% test getting all available data in the specified range from the specified
% point

  % X1 = repmat([1 2],10,1) + randn(10,2)
  X1 = [ 0.020498921332767          1.38965592083978; ...
         -0.80619687360328          2.27607461036333; ...
           1.6094931848136          1.67583675081601; ...
         0.650804938165822          2.85922989147957; ...
         -0.77928713627071          0.63905410022498; ...
          1.11199440272707          0.14438469671649; ...
         0.340288955607853          0.69203955888576; ...
         0.760824545586689          1.98345510508742; ...
           1.7390597332666          1.25858899410418; ...
           1.6610778581751          1.48258769208384];

  %{
  % X2 = repmat([-2 -1],10,1) + randn(10,2)
  X2 = [ -2.30743677363252         -2.74439891445086; ...
         -1.59642757081464         -2.07408152607031; ...
         -3.81637793040725         -1.22778729847389; ...
         -2.36725951870243        -0.390758367947946; ...
         -2.13545799887958          -1.6084521990818; ...
         -1.22393903602435          1.55523627884863; ...
         -1.20234827756917         -2.19270599748587; ...
          -1.5970690650945         -1.59555256727339; ...
         -1.08502537168358        -0.532838372665024; ...
         -1.35903245251369         -1.20420302511798];

  % X3 = repmat([-3 -3],10,1) + 2*randn(10,2)
  X3 = [ -2.71459478403865         -2.21301457334557; ...
         -2.01409121171665         -4.46756509049671; ...
         -2.06507536126567         -3.68287329420183; ...
          1.41423870282787         -3.58354233314548; ...
         0.818658559659212          -2.7094362548471; ...
         -3.71488149757023        -0.690882290497715; ...
          -3.1174898921933            -6.43533603837; ...
         -3.73171780147976         -6.95761634537453; ...
         -3.95334941789486         -5.19819827263364; ...
         -3.64829547254891         -2.92085132938035];
  %}

  xmean = [0 0];

  % f = figure();
  % scatter(X1(:,1), X1(:,2), 'r+');
  % hold on;
  % scatter(X2(:,1), X2(:,2), 'g+');
  % scatter(X3(:,1), X3(:,2), 'b+');
  % plot(xmean(1), xmean(2), 'k+');

  a = Archive(2);

  a = a.save(X1, fsphere(X1), 1);
  rangeSigma = 1;
  sigma = 1;
  BD = [1 0; 0 1.5];

  [Xt, ~] = a.getDataNearPoint(5, xmean, rangeSigma, sigma, BD);
  Xt_test = [ ...
     -0.77929      0.63905; ...
     0.020499       1.3897; ...
      0.34028      0.69204];
  [~, idx] = sort(Xt(:,1),1);
  verifyEqual(testCase, Xt(idx,:), Xt_test, 'AbsTol', 1e-2);

  % scatter(Xt(:,1), Xt(:,2), 'ro');

  % a = a.save(X2, fsphere(X2), 2);
  % sigma = 1.5;
  % [Xt, ~] = a.getDataNearPoint(5, xmean, sigma, BD);
  % [~, idx] = sort(Xt(:,1),1);
  % scatter(Xt(:,1), Xt(:,2), 'go');

end

function y = fsphere(x)
  y = sum(x.^2, 2);
end
