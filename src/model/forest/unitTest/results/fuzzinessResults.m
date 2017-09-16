load('TreeModelTestFuzziness.mat');

for fNum = [1, 2, 6, 8, 13, 14, 15, 17, 20, 21]
  name = sprintf('%02d', fNum);
  group = [];
  for i = 1:numel(results)
    if ~strcmpi(results(i).name(7:8), name)
      continue
    end
    %if isfield(results(i).params, 'modelSpec')
    if ~isfield(results(i).params, 'modelSpec') || ~strcmpi(results(i).params.modelSpec, 'quadratic')
      continue
    end
    group = [group; results(i)];
  end
  
  figure;
  params = [group.params];
  x = [params.fuzziness];
  y = [params.lambda];
  z = [group.trainRMSE];
  [~,m] = min(z);
  scatter3(x, y, z, 'b');
  hold on;
  scatter3(x(m), y(m), z(m), 'bx');
  z = [group.testRMSE];
  [~,m] = min(z);
  hold on;
  scatter3(x, y, z, 'r');
  hold on;
  scatter3(x(m), y(m), z(m), 'rx');
  title(sprintf('RMSE'));
  close;
  %figure;
  %title(sprintf('test RMSE'));
end
