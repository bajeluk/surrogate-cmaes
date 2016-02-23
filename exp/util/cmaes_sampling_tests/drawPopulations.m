function drawPopulations(FUNC, cmaOutput, cmOptions)
  
  % handle function pointers and strings with name of the function
  % correctly
  if (isstr(FUNC))
    func = str2func(FUNC);
  elseif (isa(FUNC,'function_handle'))
    func = FUNC;
  else
    error('FUNC not recognized, it must be a function_handle or a string with the function name.');
  end

  % prepare figure and axes for drawing
  f = figure();
  ax = axes();
  set(ax, 'XLim', [cmOptions.LBounds cmOptions.UBounds]);
  set(ax, 'YLim', [cmOptions.LBounds cmOptions.UBounds]);

  hold on;

  % draw contours of the fitness
  x = linspace(cmOptions.LBounds, cmOptions.UBounds, 50);
  y = linspace(cmOptions.LBounds, cmOptions.UBounds, 50);
  [X, Y] = meshgrid(x, y);
  X_lin = reshape(X, numel(X), 1);
  Y_lin = reshape(Y, numel(Y), 1);
  Z_lin = func([X_lin'; Y_lin']);
  Z = reshape(Z_lin, sqrt(length(Z_lin)), sqrt(length(Z_lin)));
  contour(X, Y, Z, 100);

  % iteratively draw the populations of each generation
  for i = 1:max(cmaOutput.generations)
    thisIdxs = cmaOutput.generations == i;
    
    % draw current population
    scatter(ax, cmaOutput.arxvalids(1,thisIdxs), cmaOutput.arxvalids(2,thisIdxs), 'bo');
    
    % identify current CMA-ES internal variables
    xmean = cmaOutput.means(:,i);   % mean 'm'
    BD = cmaOutput.BDs{i} * cmaOutput.sigmas(i);  % eig decompos. of 'C'
    
    % plot the mean 'm'
    plot(ax, xmean(1), xmean(2), 'k*');
   
    % plot the ellipsis according to 'BD'
    ellipse = plotcov2(xmean, BD*BD', 0.7);    

    title(ax, ['generation = ' num2str(i)]);
    pause;
    
    % convert the population to green
    scatter(ax, cmaOutput.arxvalids(1,thisIdxs), cmaOutput.arxvalids(2,thisIdxs), 'go');
    plot(ax, xmean(1), xmean(2), 'g*');
    
    % delete the ellipsis and the basis vectors
    delete(ellipse);
  end
  
end
