function drawPopulations(FUNC, popsCell, cmOptions)
  
  % handle function pointers and strings with name of the function
  % correctly
  if (isstr(FUNC))
    func = str2func(FUNC);
  elseif (isa(FUNC,'function_handle'))
    func = FUNC;
  else
    error('FUNC not recognized, it must be a handle or a string.');
  end

  f = figure();
  ax = axes();
  ax.XLim = [cmOptions.LBounds cmOptions.UBounds];
  ax.YLim = [cmOptions.LBounds cmOptions.UBounds];

  hold on;
  
  x = linspace(ax.XLim(1), ax.XLim(2), 50);
  y = linspace(ax.YLim(1), ax.YLim(2), 50);
  [X, Y] = meshgrid(x, y);
  X_lin = reshape(X, prod(size(X)), 1);
  Y_lin = reshape(Y, prod(size(Y)), 1);
  Z_lin = func([X_lin'; Y_lin']);
  Z = reshape(Z_lin, sqrt(length(Z_lin)), sqrt(length(Z_lin)));
  
  contour(X, Y, Z, 100);

  for i = 1:length(popsCell)
    scatter(ax, popsCell{i}(1,:), popsCell{i}(2,:), 'bo');
    pause;
    scatter(ax, popsCell{i}(1,:), popsCell{i}(2,:), 'go');
  end
  
end