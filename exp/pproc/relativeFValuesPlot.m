function handle = relativeFValuesPlot(data, varargin)
% handle = relativeFValuesPlot(data, settings)
% Plots dependences of minimal function values on function 
% evaluations / dimension for individual functions.
%
% Input:
%   data      - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataNames'     - cell array of data names (e.g. names of algorithms)
%     'DataDims'      - dimensions of data
%     'DataFuns'      - functions of data
%     'PlotDims'      - dimensions chosen to plot
%     'PlotFuns'      - functions chosen to plot
%     'Colors'        - colors of individual algorithms
%     'AggregateDims' - aggregate dimensions in plots | boolean
%     'Statistic'     - statistic of data | string or handle (@mean, 
%                       @median)
%     'MinValue'      - minimal possible function value
%     'OneFigure'     - plot in one figure | boolean
%     'SplitLegend'   - legend splitted in first two graphs | boolean 
%
% Output:
%   handle - handles of resulting figures
%
% See Also:
%   speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  if nargin < 1 || isempty(data)
    handle = [];
    help relativeFValuesPlot
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
  % TODO: rewrite to settings structure which can be given as an argument
  % to other functions
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
  aggFuns = defopts(settings, 'AggregateFuns', false);
  minValue = defopts(settings, 'MinValue', 1e-8);
  maxEval = defopts(settings, 'MaxEval', 100);
  statistic = defopts(settings, 'Statistic', @mean);
  oneFigure = defopts(settings, 'OneFigure', false);
  splitLegend = defopts(settings, 'SplitLegend', false);
  if ischar(statistic)
    if strcmp(statistic, 'quantile')
      statistic = @(x, dim) quantile(x, [0.25, 0.5, 0.75], dim);
    else
      statistic = str2func(statistic);
    end
  end

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
        data_stats{D}{f,d}(data_stats{D}{f,d} < minValue) = minValue;
      end
    end
  end
  
  % draw plot
  handle = relativePlot(data_stats, dims, BBfunc, datanames, colors, aggDims, aggFuns, maxEval, oneFigure, splitLegend);

end

function handle = relativePlot(data_stats, dims, BBfunc, datanames, colors, aggDims, aggFuns, maxEval, oneFigure, splitLegend)
% Plots scaled graph of different algorithms in all defined functions and 
% one all dimensions

  numOfData = length(data_stats);
  numOfFuncIds = length(BBfunc);
  evaldim = 1:min(length(data_stats{1}{1}), maxEval);
  minGraph = -8;
  maxGraph =  0;
  
  for f = 1:numOfFuncIds
    % find useful data and plot 
    for d = 1:length(dims)
      % find available data
      notEmptyData = true(1, numOfData);
      for dat = 1:numOfData
        notEmptyData(dat) = ~isempty(data_stats{dat}{f,d});
        if ~notEmptyData(dat)
          warning('%s is missing in function %d and dimension %d.', datanames{dat}, BBfunc(f), dims(d))
          relativeData{dat}{f, d} = [];
        end
      end
      if any(notEmptyData)
        nEmptyId = inverseIndex(notEmptyData);
        nUsefulData = sum(notEmptyData);
        % count f-values ratio
        actualData = cell2mat(arrayfun(@(D) data_stats{nEmptyId(D)}{f,d}, 1:nUsefulData, 'UniformOutput', false));
        nData = min(maxEval, size(actualData, 1));
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
  
  if aggDims
    nDimsToPlot = 1;
    for D = 1:numOfData
      for f = 1:numOfFuncIds
        nDims = size(relativeData{D}, 2);
        relativeData{D}{f, 1} = mean(cell2mat(arrayfun(@(x) relativeData{D}{f,x}, 1:nDims, 'UniformOutput', false)'));
      end
      relativeData{D} = relativeData{D}(:, 1);
    end
  else
    nDimsToPlot = length(dims);
  end
  
  if aggFuns
    nFunsToPlot = 1;
    for D = 1:numOfData
      for d = 1:length(dims)
        nFuns = size(relativeData{D}, 1);
        relativeData{D}{1, d} = mean(cell2mat(arrayfun(@(x) relativeData{D}{x,d}, 1:nFuns, 'UniformOutput', false)'), 1);
      end
      relativeData{D} = relativeData{D}(1, :);
    end
  else
    nFunsToPlot = length(BBfunc);
  end
  
  if oneFigure && (nFunsToPlot*nDimsToPlot > 1)
    nRows = ceil(nFunsToPlot*nDimsToPlot/2);
    handle = figure('Units', 'centimeters', 'Position', [1 1 20 6*nRows]);
    subplot(nRows, 2, 1)
    onePlot(relativeData, 1, 1, evaldim, maxEval, colors, ...
            datanames, aggFuns, aggDims, BBfunc, dims, true, ...
            1*splitLegend, false);
    if nDimsToPlot > 1
      f = 1;
      d = 2;
    else
      f = 2;
      d = 1;
    end
    subplot(nRows, 2, 2)
    onePlot(relativeData, f, d, evaldim, maxEval, colors, ...
            datanames, aggFuns, aggDims, BBfunc, dims, splitLegend, ...
            2*splitLegend, true);
    
    if (nFunsToPlot*nDimsToPlot > 2)
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
          omitYLabel = ~logical(mod(plotId, 2));
          subplot(nRows, 2, plotId)
          onePlot(relativeData, f, d, evaldim, maxEval, colors, ...
            datanames, aggFuns, aggDims, BBfunc, dims, false, ...
            0, omitYLabel);
        end
      end
    end
  else
    handle = zeros(1, nDimsToPlot*nFunsToPlot);
    for f = 1:nFunsToPlot
      for d = 1:nDimsToPlot
        handle((d-1) * nFunsToPlot + f) = figure('Units', 'centimeters', 'Position', [1 1 12.5 6]);
        onePlot(relativeData, f, d, evaldim, maxEval, colors, ...
                datanames, aggFuns, aggDims, BBfunc, dims, ~splitLegend || (f == 1 && d == 1), ...
                0, false);
      end
    end
  end
  
end

function onePlot(relativeData, fId, dId, evaldim, maxEval, colors, ...
                 datanames, aggFuns, aggDims, BBfunc, dims, dispLegend, ...
                 splitLegendOption, omitYLabel)
% Plots one scaled graph 
%
% Note: Omitting y-label is currently unabled. To change this status
% uncomment rows at the end of onePlot function.

  medianLineWidth = 2;

  for dat = 1:length(relativeData)
    notEmptyData(dat) = ~isempty(relativeData{dat}{fId, dId});
  end
  if any(notEmptyData)
    nEmptyId = inverseIndex(notEmptyData);
    nUsefulData = sum(notEmptyData);
    h = zeros(1, nUsefulData);
    ftitle = cell(1, nUsefulData);
    % mean
    h(1) = plot(evaldim, relativeData{nEmptyId(1)}{fId, dId}(1:maxEval), ...
      'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(1), :));
    ftitle{1} = datanames{nEmptyId(1)};
    hold on
    grid on
    for dat = 2:nUsefulData
      h(dat) = plot(evaldim, relativeData{nEmptyId(dat)}{fId, dId}(1:maxEval), ...
        'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(dat), :));
      ftitle{dat} = datanames{nEmptyId(dat)};
    end
    if dispLegend
      switch splitLegendOption
        case 0
          legIds = true(1, nUsefulData);
        case 1
          legIds = [true(1, floor(nUsefulData/2)), false(1, nUsefulData - floor(nUsefulData/2))];
        case 2
          legIds = [false(1, floor(nUsefulData/2)), true(1, nUsefulData - floor(nUsefulData/2))];
      end
      legend(h(legIds), ftitle(legIds), 'Location', 'NorthEast')
    end
  else
    warning('Function %d dimension %d has no data available', BBfunc(fId), dims(dId))
  end

  titleString = '';
  if ~aggFuns
    titleString = ['f', num2str(BBfunc(fId))];
  end
  if ~aggDims
    titleString = [titleString, ' ', num2str(dims(dId)),'D'];
  end
  title(titleString)
  xlabel('Number of evaluations / D')
%   if ~omitYLabel
    ylabel('\Delta_f^{log}')
%   end
  hold off

end
