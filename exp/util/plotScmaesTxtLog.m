%PLOTSCMAESTXTLOG Plots variables from S-CMA-ES log in txt logfile
%
% txt logfiles are located in bbob_output/*_log_*.dat
%
% Optional parameter--value settings:
% 'xAxis'       'iters'  -- CMA-ES iterations on x-axis
%               'fevals' -- original function evaluations on x-axis
%
% 'instance'    1,2,...  -- plot specified instance (of this order) instead 
%                           of the first instance in the txt log file
%
% 'title'       [string] -- string which shall be placed into title instaed
%                           of the default one (f, D, EXPID, ID, instance)
function fig = plotScmaesTxtLog(exp_id, fun, dim, id, varargin)
  % parse optional input arguments
  if (~isempty(varargin))
    if (~isstruct(varargin))
      try
        opts = struct(varargin{:});
      catch
      end
    else
      opts = varargin{1};
    end
  else
    opts = struct();
  end
  opts.xAxis = defopts(opts, 'xAxis', 'iters');
  opts.instance = defopts(opts, 'instance', 1);

  % plot the plots for each instance
  for inst = opts.instance

  thisTitle = defopts(opts, 'title', sprintf('f%d in %dD  %s  ID=%02d  inst=%d', ...
    fun, dim, strrep(exp_id, '_', '-'), id, inst));

  % load data from txt log file
  t = readScmaesTxtLog(exp_id, fun, dim, id, inst);

  fig = figure();
  fig.Position = [50, 50, 1200, 550];

  % determine the x-axis
  if (strcmp(opts.xAxis, 'iters'))
    xVals = correctItersAfterRestarts(t.iter);
    xLabel = 'CMA-ES iterations';
  elseif (strcmp(opts.xAxis, 'fevals'))
    xVals = t.totEvals;
    xLabel = 'function evaluations';
  else
    error(sprintf('X Axis type is not supported: %s', opts.xAxis));
  end

  % prepare plot legend strings
  legends = {};

  % plot the graphs corresponding to the left y-axis
  semilogy(xVals, t.sigmaPow2, 'Color', 'blue', 'LineWidth', 2);
  legends{end+1} = 'sigma^2';
  hold on;
  semilogy(xVals, t.Dopt, 'Color', 'black', 'LineWidth', 2.5);
  legends{end+1} = 'Dopt';
  semilogy(xVals, t.rmseReeval)
  legends{end+1} = 'RMSE on re-evaled';
  semilogy(xVals, t.diagD_ratio);
  legends{end+1} = 'max/min of diag(D)';

  % plot graphs of diagonal elements of CMA-ES covariance matrix
  fieldnames = cell(1,dim);
  for d = 1:dim
    fieldnames{d} = sprintf('diagD_%02d', d);
    semilogy(xVals, table2array(t(:,fieldnames{d})));
    legends{end+1} = strrep(fieldnames{d}, '_', '-');;
  end

  % switch to the right y-axis
  yyaxis right

  plot(xVals, t.nKendallOldModel, 'r--');
  legends{end+1} = 'Kendall of old model';
  plot(xVals, t.rnkValid, 'Color', [197 217 66]./255);  % yellow-green
  legends{end+1} = 'rankDiff on validation set';

  % adaptive statistics
  if (any(strcmp('adaptSmoothedErr', t.Properties.VariableNames)))
    plot(xVals, t.adaptSmoothedErr, 'g-', 'LineWidth', 1.5);
    legends{end+1} = 'Smoothed Err';
  end
  if (any(strcmp('adaptErr', t.Properties.VariableNames)))
    plot(xVals, t.adaptErr, 'Color', [0 150 0]./255, 'LineStyle', '-'); % dark green
    legends{end+1} = 'Actual Err';
  end

  % original FE's in generations
  oFE = t.orEvalsPop ./ max(t.orEvalsPop);
  plot(xVals, oFE, 'm-', 'LineWidth', 1.5);
  legends{end+1} = '# orig. evals';

  % fill the legend and axis labels
  legend(legends{:});
  xlabel(xLabel);
  ylabel(sprintf('# orig. evals: 1.0 == %d points', max(t.orEvalsPop)));
  title(thisTitle);

  end % for inst = opts.instance

end

% correct iteration numbers after restarts (it goes wrongly down to zero)
function corrIters = correctItersAfterRestarts(iters)
  corrIters = zeros(size(iters));
  last = 0;
  fix  = 0;
  for i = 1:length(iters)
    if (iters(i) < last)
      fix = fix + last;
    end
    corrIters(i) = fix + iters(i);
    last = iters(i);
  end
end
