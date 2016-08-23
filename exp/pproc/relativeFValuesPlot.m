function handle = relativeFValuesPlot(data, varargin)
% handle = relativeFValuesPlot(data, settings)
% Plots dependences of minimal function values on function 
% evaluations / dimension for individual functions.
%
% Input:
%   data     - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'AggregateDims' - aggregate dimensions in plots | boolean
%     'AggregateFuns' - aggregate functions in plots | boolean
%     'Colors'        - vector Nx3 of (RGB) colors of individual 
%                       algorithms, N is the number of algorithms | double 
%                       array
%                       (e.g. 2 algorithms: [115, 60, 215; 0, 40, 60]) 
%     'DataNames'     - cell array of data names (e.g. names of algorithms)
%                       | cell array of strings
%                       (e.g. {'alg1', 'alg2', 'alg3'})
%     'DataDims'      - dimensions of data | integer vector 
%                       (e.g. [2, 3, 5, 10])
%     'DataFuns'      - functions of data | integer vector
%                       (e.g. [1, 2, 3, 4, 5, 11, 12])
%     'LegendOption'  - legend settings:
%                         'show'  - show legend
%                         'hide'  - do not show legend
%                         'split' - legend splitted in first two graphs
%                         'out'   - legend is out of the last graph
%     'LineSpec'      - specification of plotted lines (see help plot ->
%                       LineSpec), to set colors use 'Colors' settings | 
%                       cell array of string
%                       (e.g. 3 algorithms: {'-.', '-o', '.'})
%     'MaxEval'       - maximum of evaluations divided by dimension |
%                       integer
%     'MinValue'      - minimal possible function value | double
%     'OmitYLabel'    - omit y-label in even plots | boolean
%     'OneFigure'     - all plots in one figure | boolean
%     'PlotDims'      - dimensions chosen to plot | integer vector
%                       (e.g. [5, 10])
%     'PlotFuns'      - functions chosen to plot | integer vector
%                       (e.g. [3, 6, 10, 17])
%     'Statistic'     - statistic of data | string or handle (@mean, 
%                       @median)
%
% Output:
%   handle - handles of resulting figures
%
% See Also:
%   speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  if nargin < 1 || isempty(data)
    if nargout > 0
      handle = [];
    end
    help relativeFValuesPlot
    return
  end
  settings = settings2struct(varargin);
  
  % parse settings
  numOfData = length(data);
  plotSet.datanames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:numOfData, 'UniformOutput', false));
  assert(length(plotSet.datanames) == numOfData, 'Number of data and number of DataNames are not the same')
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  plotSet.dims    = defopts(settings, 'PlotDims', funcSet.dims);
  plotSet.BBfunc  = defopts(settings, 'PlotFuns', funcSet.BBfunc);
  colors  = defopts(settings, 'Colors', rand(numOfData, 3));
  if max(colors) > 1
    colors = colors / 255;
  end
  plotSet.colors = colors;
  plotSet.aggDims = defopts(settings, 'AggregateDims', false);
  plotSet.aggFuns = defopts(settings, 'AggregateFuns', false);
  plotSet.minValue = defopts(settings, 'MinValue', 1e-8);
  plotSet.maxEval = defopts(settings, 'MaxEval', 250);
  statistic = defopts(settings, 'Statistic', @mean);
  plotSet.oneFigure = defopts(settings, 'OneFigure', false);
  plotSet.legendOption = defopts(settings, 'LegendOption', 'show');
  defaultLine = arrayfun(@(x) '-', 1:numOfData, 'UniformOutput', false);
  plotSet.lineSpec = defopts(settings, 'LineSpec', defaultLine);
  if length(plotSet.lineSpec) ~= numOfData
    warning('Number of line specification strings and number of data are different. Setting default values.')
    plotSet.lineSpec = defaultLine;
  end
  plotSet.omitYLabel = defopts(settings, 'OmitYLabel', false);
  if ischar(statistic)
    if strcmp(statistic, 'quantile')
      statistic = @(x, dim) quantile(x, [0.25, 0.5, 0.75], dim);
    else
      statistic = str2func(statistic);
    end
  end

  % get function and dimension IDs
  dimInvIds = inverseIndex(funcSet.dims);
  dimIds = dimInvIds(plotSet.dims);
  funcInvIds = inverseIndex(funcSet.BBfunc);
  funcIds = funcInvIds(plotSet.BBfunc);

  if ~all(dimIds)
    fprintf('Wrong dimesion request!\n')
  end
  if ~all(funcIds)
    fprintf('Wrong function request!\n')
  end

  % count means
  useMaxInstances = 15;
  data_stats = cellfun(@(D) gainStatistic(D, dimIds, funcIds, ...
                            'MaxInstances', useMaxInstances, ...
                            'AverageDims', false, ...
                            'Statistic', statistic, ...
                            'SuppWarning', true), ...
                       data, 'UniformOutput', false);
                          
  % minimal value cannot be lower than minValue
  for D = 1:numOfData
    for f = 1:length(funcIds)
      for d = 1:length(dimIds)
        data_stats{D}{f,d}(data_stats{D}{f,d} < plotSet.minValue) = plotSet.minValue;
      end
    end
  end
  
  % draw plot
  handle = relativePlot(data_stats, plotSet);

end

function handle = relativePlot(data_stats, settings)
% Plots scaled graph of different algorithms in all defined functions and 
% one all dimensions

  numOfData = length(data_stats);
  numOfFuncIds = length(settings.BBfunc);
  minGraph = -8;
  maxGraph =  0;
  dispLegend = true;
  splitLegend = false;
  if strcmp(settings.legendOption, 'split')
    splitLegend = true;
  end
  if strcmp(settings.legendOption, 'hide')
    dispLegend = false;
  end
  settings.legendLocation = 'NorthEast';

  
  for f = 1:numOfFuncIds
    % find useful data and plot 
    for d = 1:length(settings.dims)
      % find available data
      notEmptyData = true(1, numOfData);
      for dat = 1:numOfData
        notEmptyData(dat) = ~isempty(data_stats{dat}{f,d});
        if ~notEmptyData(dat)
          warning('%s is missing in function %d and dimension %d.', settings.datanames{dat}, settings.BBfunc(f), settings.dims(d))
          relativeData{dat}{f, d} = [];
        end
      end
      if any(notEmptyData)
        nEmptyId = inverseIndex(notEmptyData);
        nUsefulData = sum(notEmptyData);
        % count f-values ratio
        actualData = cell2mat(arrayfun(@(D) data_stats{nEmptyId(D)}{f,d}, 1:nUsefulData, 'UniformOutput', false));
        nData = min(settings.maxEval, size(actualData, 1));
        actualData = log10( actualData(1:nData, :) );
        actualMin = min(min(actualData));
        actualMax = max(max(actualData));
        for D = 1:nUsefulData
          % % this is old version for scaling BEFORE logarithm
          % minGraph = 1e-8; maxGraph = 1;
          % relativeData{nEmptyId(D)}{f, d} = log10(thisData);
          thisData = ((actualData(:, D) - actualMin) * (maxGraph - minGraph)/(actualMax - actualMin))' + minGraph;
          relativeData{nEmptyId(D)}{f, d} = thisData;
        end
      end
    end
  end
  
  % aggregate accross dimensions
  if settings.aggDims
    nDimsToPlot = 1;
    for D = 1:numOfData
      for f = 1:numOfFuncIds
        nDims = size(relativeData{D}, 2);
        relativeData{D}{f, 1} = mean(cell2mat(arrayfun(@(x) relativeData{D}{f,x}, 1:nDims, 'UniformOutput', false)'));
      end
      relativeData{D} = relativeData{D}(:, 1);
    end
  else
    nDimsToPlot = length(settings.dims);
  end
  
  % aggregate accross functions
  if settings.aggFuns
    nFunsToPlot = 1;
    for D = 1:numOfData
      for d = 1:length(settings.dims)
        nFuns = size(relativeData{D}, 1);
        relativeData{D}{1, d} = mean(cell2mat(arrayfun(@(x) relativeData{D}{x,d}, 1:nFuns, 'UniformOutput', false)'), 1);
      end
      relativeData{D} = relativeData{D}(1, :);
    end
  else
    nFunsToPlot = length(settings.BBfunc);
  end
  
  nPlots = nFunsToPlot*nDimsToPlot;
  
  % all plots one figure
  if settings.oneFigure && (nPlots > 1)
    nRows = ceil(nPlots/2);
    handle = figure('Units', 'centimeters', 'Position', [1 1 20 6*nRows]);
    % first plot
    subplot(nRows, 2, 1)
    onePlot(relativeData, 1, 1, settings, dispLegend, ...
            1*splitLegend, false);
    if nDimsToPlot > 1
      f = 1;
      d = 2;
    else
      f = 2;
      d = 1;
    end
    % second plot
    subplot(nRows, 2, 2)
    onePlot(relativeData, f, d, settings, splitLegend, ...
            2*splitLegend, settings.omitYLabel);
    % the rest of plots
    if (nPlots > 2)
      if nDimsToPlot > 2
        fStart = 1;
        dStart = 3;
      else
        fStart = f + 1;
        dStart = 1;
      end

      for f = fStart:nFunsToPlot
        for d = dStart:nDimsToPlot
          plotId = (d-1) * nFunsToPlot + f;
          omitYLabelStatus = ~logical(mod(plotId, 2));
          subplot(nRows, 2, plotId)
          onePlot(relativeData, f, d, settings, false, ...
            0, settings.omitYLabel && omitYLabelStatus);
        end
      end
    end
    
  % one plot one figure
  else
    
    % legend settings
    if strcmp(settings.legendOption, 'out')
      handle = cell(1, nPlots + 1);
      legendFigNum = nPlots + 1;
      settings.legendLocation = 'EastOutside';
    else
      handle = cell(1, nPlots);
      legendFigNum = 1;
    end
    
    % plot all functions and dimensions
    plottedInAny = false(1, numOfData);
    for f = 1:nFunsToPlot
      for d = 1:nDimsToPlot
        handle{(d-1) * nFunsToPlot + f} = ...
          figure('Units', 'centimeters', 'Position', [1, 1, 12.5, 6]);
        plottedInAny = plottedInAny | ...
          onePlot(relativeData, f, d, settings, dispLegend && (strcmp(settings.legendOption, 'show') || (f*d == legendFigNum)), ...
                0, false);
      end
    end
    if strcmp(settings.legendOption, 'out')
      handle{legendFigNum} = soloLegend(settings.colors(plottedInAny, :), settings.datanames(plottedInAny), 2);
    end
  end
  
end

function notEmptyData = onePlot(relativeData, fId, dId, ...
                 settings, dispLegend, splitLegendStatus, omitYLabel)
% Plots one scaled graph 
%
% Note: Omitting y-label is currently enabled. To change this status
% uncomment rows at the end of onePlot function.

  % parsing settings
  maxEval = settings.maxEval;
  colors = settings.colors;
  datanames = settings.datanames;
  aggFuns = settings.aggFuns;
  aggDims = settings.aggDims;
  BBfunc = settings.BBfunc;
  dims = settings.dims;
  lineSpec = settings.lineSpec;

  % default value
  medianLineWidth = 1;

  % find not empty data
  for dat = 1:length(relativeData)
    notEmptyData(dat) = ~isempty(relativeData{dat}{fId, dId});
  end
  
  if any(notEmptyData)
    nEmptyId = inverseIndex(notEmptyData);
    nUsefulData = sum(notEmptyData);
    % find minimal number of function evaluations to plot
    evaldim = 1:min([arrayfun(@(x) length(relativeData{nEmptyId(x)}{fId, dId}), 1:nUsefulData), maxEval]);
    h = zeros(1, nUsefulData);
    ftitle = cell(1, nUsefulData);
    % plot first line 
    h(1) = plot(evaldim, relativeData{nEmptyId(1)}{fId, dId}(evaldim), ...
      lineSpec{1}, 'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(1), :));
    ftitle{1} = datanames{nEmptyId(1)};
    hold on
    grid on
    % plot rest of lines
    for dat = 2:nUsefulData
      h(dat) = plot(evaldim, relativeData{nEmptyId(dat)}{fId, dId}(evaldim), ...
        lineSpec{dat}, 'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(dat), :));
      ftitle{dat} = datanames{nEmptyId(dat)};
    end
    % legend settings
    if dispLegend
      switch splitLegendStatus
        case 0
          legIds = true(1, nUsefulData);
        case 1
          legIds = [true(1, floor(nUsefulData/2)), false(1, nUsefulData - floor(nUsefulData/2))];
        case 2
          legIds = [false(1, floor(nUsefulData/2)), true(1, nUsefulData - floor(nUsefulData/2))];
      end
      legend(h(legIds), ftitle(legIds), 'Location', settings.legendLocation)
    end
  else
    warning('Function %d dimension %d has no data available', BBfunc(fId), dims(dId))
  end

  % make title
  titleString = '';
  if ~aggFuns
    titleString = ['f', num2str(BBfunc(fId))];
  end
  if ~aggDims
    titleString = [titleString, ' ', num2str(dims(dId)),'D'];
  end
  title(titleString)
  xlabel('Number of evaluations / D')
  if ~omitYLabel
    ylabel('\Delta_f^{log}')
  end
  
  hold off

end

function h = soloLegend(colors, names, lineWidth)
% plots legend in solo figure

  nNames = length(names);
  yDiff = 1/(2*nNames);
  margin = yDiff;
  x_line_length = 2*margin;
  
  h = figure('Units', 'centimeters', 'Position', [1, 1, 12.5, 6]);
  hold on
  axis off
  % first frame point
  fa = [0, 0]; % [fa_x, fa_y]
  
  % lines and text
  x_line_coor = fa(1) + margin + x_line_length;
  font = min(10, 125/nNames);
  for n = 1:nNames
    y_coor = -2*n*yDiff;
    line([fa(1) + margin, x_line_coor], y_coor*[1,1], ...
         'Color', colors(n, :), 'LineWidth', lineWidth)
    text(x_line_coor + margin, y_coor, names{n}, 'FontSize', font)
  end
  
  % second frame point
  fb = [1, -2*(nNames+1)*yDiff]; % [fb_x, fb_y]
  line([fa(1), fb(1)], [fa(2), fa(2)], 'Color', 'k') %  _
  line([fa(1), fa(1)], [fa(2), fb(2)], 'Color', 'k') % |
  line([fb(1), fb(1)], [fa(2), fb(2)], 'Color', 'k') %  _|
  line([fa(1), fb(1)], [fb(2), fb(2)], 'Color', 'k') %  
  
  hold off
end