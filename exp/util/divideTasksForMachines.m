function [c, usedSizes] = divideTasksForMachines(nMachines, dimensions, func)
  % solves an n-subset problem: divide a vector v = func(dimensions)
  % into nMachines sets in order to have the same sums
  c = cell(1,nMachines);
  sizes = func(dimensions);
  usedSizes = zeros(1,nMachines);

  [sortSizes, sizesIdxs] = sort(sizes, 2, 'descend');
  for i = 1:length(sizesIdxs)
    [~, minMachine] = min(usedSizes);
    c{minMachine} = [sizesIdxs(i) c{minMachine}];
    usedSizes(minMachine) = usedSizes(minMachine) + sortSizes(i);
  end
end

