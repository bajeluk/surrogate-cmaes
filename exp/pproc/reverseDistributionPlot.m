function handle = reverseDistributionPlot(data, varargin)
% handle = reverseDistributionPlot(data, settings)

  % initialization
  if nargin < 1 || isempty(data)
    handle = [];
    help reverseDistributionPlot
    return
  end
  if isstruct(varargin)
    settings = varargin;
  else
    % keep cells as cells due to struct command
    vararCellId = cellfun(@iscell, varargin);
    varargin(vararCellId) = {varargin(vararCellId)};
    settings = struct(varargin{:});
  end
  
  % initalize settings
  numOfData = length(data);
  datanames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:numOfData, 'UniformOutput', false));
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  dims    = defopts(settings, 'PlotDims', funcSet.dims);
  BBfunc  = defopts(settings, 'PlotFuns', funcSet.BBfunc);
  colors  = defopts(settings, 'Colors', rand(numOfData, 3));
  aggDims = defopts(settings, 'AggregateDims', false);
  defTargets = defopts(settings, 'DefaultTargets', [2*(1:25), 5*(11:20), 10*(11:25)]);
  useMaxEvals = defopts(settings, 'MaxEvals', 250);
  
  % get function and dimension IDs
  dimInvIds = inverseIndex(funcSet.dims);
  dimIds = dimInvIds(dims);
  funcInvIds = inverseIndex(funcSet.BBfunc);
  funcIds = funcInvIds(BBfunc);
  
  if ~all(dimIds)
    fprintf('Wrong dimesion request!\n')
  end
  if ~all(funcIds)
    fprintf('Wrong function request!\n')
  end
  
  % count relative distances
  useMaxInstances = 15;
  [targetEvals, yTargets] = ...
    relativeMeasure(data, dimIds, funcIds, ...
                   'DefaultTargets', defTargets, 'AggregateDims', aggDims, ...
                   'RefMeasure', @(x, y) min(x, [], y), 'MaxInstances', useMaxInstances);
  
  % plot values
  lineWidth = 2;
  
  if aggDims
    nDimsToPlot = 1;
  else
    nDimsToPlot = length(dims);
  end
  numOfFuncIds = length(funcIds);
  
  handle = zeros(1, nDimsToPlot*numOfFuncIds);
  for d = 1:nDimsToPlot
    for f = 1:numOfFuncIds
      handle((d-1) * numOfFuncIds + f) = figure('Units', 'centimeters', 'Position', [1 1 12.5 6]);
      % find available data
      notEmptyData = true(1, numOfData);
      for dat = 1:numOfData
        notEmptyData(dat) = ~isempty(targetEvals{dat}{f,d});
      end
      if any(notEmptyData)
        nEmptyId = inverseIndex(notEmptyData);
        h = zeros(1, sum(notEmptyData));
        ftitle = cell(1, sum(notEmptyData));
        
        % gain values to plot
        allPlotValues = cell2mat(arrayfun(@(x) (targetEvals{x}{f,d}), 1:sum(notEmptyData), 'UniformOutput', false)');
        numOfEvals = min([find(all(allPlotValues == max(max(allPlotValues))), 1, 'first'), useMaxEvals]);
        numOfTargets = length(yTargets{f,d});
        plotValues = numOfTargets - allPlotValues(1, 1:numOfEvals);
        h(1) = plot(1:numOfEvals, plotValues, ...
          'LineWidth', lineWidth, 'Color', colors(nEmptyId(1), :));
        ftitle{1} = datanames{nEmptyId(1)};
        hold on
        for dat = 2:sum(notEmptyData)
          plotValues = numOfTargets - allPlotValues(dat, 1:numOfEvals);
          h(dat) = plot(1:numOfEvals, plotValues, ...
          'LineWidth', lineWidth, 'Color', colors(nEmptyId(dat), :));
          ftitle{dat} = datanames{nEmptyId(dat)};
        end
        legend(h, ftitle, 'Location', 'NorthEast')
      else
        warning('Function %d dimension %d has no data available', BBfunc(f), dims(d))
      end
      
      % final title and axis labeling
      if aggDims
        title(['f', num2str(BBfunc(f))])
      else
        title(['f', num2str(BBfunc(f)), ' ', num2str(dims(d)),'D'])
      end
      ax = gca;
      if ~aggDims
        reverseYTargets = yTargets{f,d}(end:-1:1);
        ax.YTickLabels = [arrayfun(@(x) sprintf('%0.1e', x), reverseYTargets(ax.YTick(1:end-1) + 1)', 'UniformOutput', false), {''}];
      end
      xlabel('Number of evaluations / D')
      ylabel('Targets')
      hold off
    end
  end
end