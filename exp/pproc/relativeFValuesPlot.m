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
%     'FunctionNames' - show function names in header | boolean
%     'LegendOption'  - legend settings:
%                         'show'    - show legend
%                         'hide'    - do not show legend
%                         'first'   - ledend only in the first graph
%                         'split'   - legend splitted in first two graphs
%                         'out'     - legend is in one separate figure
%                         'manyout' - legend is in multiple separated
%                                     figures
%     'LineSpec'      - specification of plotted lines (see help plot ->
%                       LineSpec), to set colors use 'Colors' settings | 
%                       cell array of string
%                       (e.g. 3 algorithms: {'-.', '-o', '.'})
%     'LineWidth'     - line width of individual lines | double
%     'MaxEval'       - maximum of evaluations divided by dimension |
%                       integer
%     'MinValue'      - minimal possible function value | double
%     'OmitEmpty'     - omit plots with no data available| boolean
%     'OmitYLabel'    - omit y-label in even plots | boolean
%     'OneFigure'     - all plots in one figure | boolean
%     'PlotDims'      - dimensions chosen to plot | integer vector
%                       (e.g. [5, 10])
%     'PlotFuns'      - functions chosen to plot | integer vector
%                       (e.g. [3, 6, 10, 17])
%     'Statistic'     - statistic of data | string or handle (@mean, 
%                       @median)
%     'PlotGrid'      - the number of plots in the final plot grid in order
%                       to omit the y-label or x-label in internal plots
%                       2-element vector of number of plots in y and x axis,
%                       or [] for plotting x/y-labels everywhere
%     'ScaleY08'      - whether to always scale the y-axis to [-8, 0]
%                       true by default
%
% Output:
%   handle - handles of resulting figures | cell-array
%
% See Also:
%   speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  if nargin < 1 || isempty(data)
    if nargout > 0
      handle = {};
    end
    help relativeFValuesPlot
    return
  end
  settings = settings2struct(varargin);
  
  % checkout data
  emptyData = cellfun(@isempty, data);
  if any(emptyData)
    emptyDataString = num2str(find(emptyData), '%d, ');
    warning(['Data cells %s are empty and will not be plotted and data-dependent ', ...
             'properties may not be accurate, e.g. Colors, LineSpec, etc.'], ...
             emptyDataString(1:end-1))
    % TODO: exclude missing data from all properties
    data = data(~emptyData);
  elseif all(emptyData)
    error('All data cells are empty')
  end
  
  % parse settings
  numOfData = length(data);
  % name settings
  plotSet.datanames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:numOfData, 'UniformOutput', false));
  if length(plotSet.datanames) ~= numOfData && any(emptyData)
    plotSet.datanames = plotSet.datanames(~emptyData);
  end
  assert(length(plotSet.datanames) == numOfData, 'Number of data and number of DataNames are not the same')
  % function and dimension settings
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  plotSet.dims    = defopts(settings, 'PlotDims', funcSet.dims);
  plotSet.BBfunc  = defopts(settings, 'PlotFuns', funcSet.BBfunc);
  plotSet.aggDims = defopts(settings, 'AggregateDims', false);
  plotSet.aggFuns = defopts(settings, 'AggregateFuns', false);
  % color settings
  colors  = defopts(settings, 'Colors', rand(numOfData, 3));
  if max(colors) > 1
    colors = colors / 255;
  end
  if size(colors, 1) ~= numOfData && any(emptyData)
    colors = colors(~emptyData, :);
  end
  plotSet.colors = colors;
  % plot settings
  plotSet.funcNames = defopts(settings, 'FunctionNames', false);
  plotSet.minValue = defopts(settings, 'MinValue', 1e-8);
  plotSet.maxEval = defopts(settings, 'MaxEval', 250);
  plotSet.oneFigure = defopts(settings, 'OneFigure', false);
  plotSet.legendOption = lower(defopts(settings, 'LegendOption', 'show'));
  plotSet.omitEmptyPlots = defopts(settings, 'OmitEmpty', false);
  plotSet.omitYLabel = defopts(settings, 'OmitYLabel', false);
  plotSet.plotGrid = defopts(settings, 'PlotGrid', []);
  plotSet.scaleY08 = defopts(settings, 'ScaleY08', true);
  % plot-line settings
  defaultLine = arrayfun(@(x) '-', 1:numOfData, 'UniformOutput', false);
  plotSet.lineSpec = defopts(settings, 'LineSpec', defaultLine);
  plotSet.drawQuantiles = defopts(settings, 'Quantiles', false(1, numOfData));
  plotSet.TitleString = defopts(settings, 'TitleString', '');
  if length(plotSet.lineSpec) ~= numOfData && any(emptyData)
    plotSet.lineSpec = plotSet.lineSpec(~emptyData);
  end
  if length(plotSet.lineSpec) ~= numOfData
    warning('Number of line specification strings and number of data are different. Setting default values.')
    plotSet.lineSpec = defaultLine;
  end
  plotSet.lineWidth = defopts(settings, 'LineWidth', 1);
  plotSet.markers   = defopts(settings, 'Markers', {''});
  % statistic settings
  statistic = defopts(settings, 'Statistic', @mean);
  if ischar(statistic)
    if strcmp(statistic, 'quantile')
      statistic = @(x, dim) quantile(x, [0.5, 0.25, 0.75], dim);
      plotSet.lineWidth = defopts(settings, 'LineWidth', 2);
      if (size(plotSet.drawQuantiles, 2) ~= numOfData)
        warning('Quantiles can be either scalar boolean, or vector of bools with the same length as the number of data. Using the first value.');
        plotSet.drawQuantiles = plotSet.drawQuantiles(1);
      end
      if (length(plotSet.drawQuantiles) == 1)
        plotSetQuantiles = repmat(plotSet.drawQuantiles, 1, numOfData);
      end
      if (~any(plotSet.drawQuantiles))
        plotSet.drawQuantiles = true(1, numOfData);
      end
    else
      % other than quantile statistics
      if (any(plotSet.drawQuantiles))
        warning('Quantiles can be plotted only if Statistic == ''quantile''.');
        plotSet.drawQuantiles = false(1, numOfData);
      end
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
  splitLegend = false;
  if strcmp(settings.legendOption, 'split')
    splitLegend = true;
  end
  settings.legendLocation = 'NorthEast';

  if (any(settings.drawQuantiles))
    % copy the settings in order to have settings for quartiles, too
    if (size(settings.lineWidth, 2) == numOfData)
      settings.lineWidth = repmat(settings.lineWidth, 1, 3);
    else
      settings.lineWidth = repmat(settings.lineWidth, 1, 3*numOfData);
    end
    if (size(settings.lineSpec, 2) == numOfData)
      settings.lineSpec = repmat(settings.lineSpec, 1, 3);
    else
      settings.lineSpec = repmat(settings.lineSpec, 1, 3*numOfData);
    end
    settings.lineSpec((numOfData+1):end) = repmat({':'}, 1, 2*numOfData);
    if (size(settings.datanames,2) == numOfData)
      settings.datanames = repmat(settings.datanames, 1, 3);
    else
      error('The number of datanames is not consistent.');
    end
    % set up markers only for median lines, not for quantiles:
    if (length(settings.markers) ~= numOfData)
      settings.markers = repmat(settings.markers(1), 1, numOfData);
    end
    settings.markers((numOfData+1):(3*numOfData)) = repmat({''}, 1, 2*numOfData);
    
    if (size(settings.colors, 1) == numOfData)
      settings.colors = repmat(settings.colors, 3, 1);
    else
      settings.colors = repmat(settings.colors, 3*numOfData, 1);
    end
  end
  
  relativeData = cell(1, numOfData);
  % function loop
  for f = 1:numOfFuncIds
    % find useful data and plot 
    for d = 1:length(settings.dims)
      % find available data
      isNotEmptyData = ~arrayfun(@(dat) isempty(data_stats{dat}{f,d}), 1:numOfData);
      if ~all(isNotEmptyData) && any(isNotEmptyData)
        warning('%s are missing in function %d dimension %d.', ...
          strjoin(settings.datanames(~isNotEmptyData), ', '), settings.BBfunc(f), settings.dims(d))
      end
      % assign empty set to empty data
      for dat = 1:numOfData
        if ~isNotEmptyData(dat)
          relativeData{dat}{f, d} = [];
          % this data definitely wouldn't have quantile plots
          settings.drawQuantiles(dat) = false;
        end
      end
      if any(isNotEmptyData)
        if (any(settings.drawQuantiles) && length(settings.drawQuantiles) ~= length(isNotEmptyData))
          warning('''Quantiles'' parameter has to be of the same length as the number of data.');
          settings.drawQuantiles = false(1, numOfData);
        end
        nonEmptyId = find(isNotEmptyData);

        % count f-values ratio
        nData = settings.maxEval;
        
        if (any(settings.drawQuantiles))
            actualData = NaN(nData, numel(nonEmptyId)*3);
        else
            actualData = NaN(nData, numel(nonEmptyId));
        end
        % gather actual (function, dimension) non-empty data and fill the
        % missing values to maxEvals with the best achieved result
        for dat = 1:numel(nonEmptyId)
          neDat = nonEmptyId(dat);
          ds_act = data_stats{neDat}{f, d};
          if numel(ds_act) < nData
            actualData(:, dat) = [ds_act; ds_act(end)*ones(nData - numel(ds_act), 1)];
          else
              if (any(settings.drawQuantiles))
                  actualData(:, ((dat-1)*3)+1:((dat-1)*3)+3) = ds_act(1:nData, :);
              else
                  actualData(:, dat) = ds_act(1:nData);
              end
          end
        end
        % logarithmic scaling
        actualData = log10( actualData(1:nData, :) );
        actualMin = min(min(actualData));
        actualMax = max(max(actualData));

        if (any(settings.drawQuantiles))
          nonEmptyId = [nonEmptyId, numOfData+nonEmptyId, 2*numOfData+nonEmptyId];
          actualData = actualData(:,[1:3:end, 2:3:end, 3:3:end]);
          nMedians = size(actualData,2)/3;
          actualMin = min(min(actualData(:, [true(1,nMedians), repmat(settings.drawQuantiles, 1, 2)])));
          actualMax = max(max(actualData(:, [true(1,nMedians), repmat(settings.drawQuantiles, 1, 2)])));
        end

        for D = 1:size(actualData, 2)
          % % this is old version for scaling BEFORE logarithm
          % minGraph = 1e-8; maxGraph = 1;
          % relativeData{nEmptyId(D)}{f, d} = log10(thisData);
          if (settings.scaleY08 || (settings.aggDims || settings.aggFuns))
            thisData = ((actualData(:, D) - actualMin) * (maxGraph - minGraph)/(actualMax - actualMin))' + minGraph;
          else
            thisData = actualData(:, D);
          end
          relativeData{nonEmptyId(D)}{f, d} = thisData;
        end
      end
    end
  end
  
  % aggregate accross dimensions
  if settings.aggDims
    nDimsToPlot = 1;
    for D = 1:length(relativeData)
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
    actualDisp = dispLegend(1, settings.legendOption);
    onePlot(relativeData, 1, 1, settings, actualDisp, ...
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
    actualDisp = dispLegend(2, settings.legendOption);
    onePlot(relativeData, f, d, settings, actualDisp, ...
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
          plotId = (f-1) * nDimsToPlot + d;
          omitYLabelStatus = ~logical(mod(plotId, 2));
          subplot(nRows, 2, plotId)
          actualDisp = dispLegend(plotId, settings.legendOption);
          onePlot(relativeData, f, d, settings, actualDisp, ...
            0, settings.omitYLabel && omitYLabelStatus);
        end
      end
    end
    
  % one plot one figure
  else
    
    % legend settings
    if any(strcmp(settings.legendOption, {'out', 'manyout'}))
      handle = cell(1, nPlots + 1);
      settings.legendLocation = 'EastOutside';
    else
      handle = cell(1, nPlots);
    end
    
    % plot all functions and dimensions
    plottedInAny = false(1, numOfData);
    handleId = 1;
    for f = 1:nFunsToPlot
      for d = 1:nDimsToPlot
        % increase handleId if empty plots are not omitted
        if ~settings.omitEmptyPlots
          handleId = (f-1) * nDimsToPlot + d;
        end
        handle{handleId} = ...
          figure('Units', 'centimeters', 'Position', [1, 1, 12.5, 6]);
        % display legend indicator 
        actualDisp = dispLegend(handleId, settings.legendOption);
        if (~isempty(settings.plotGrid) && length(settings.plotGrid) == 2)
          omitXLabel = (mod(floor((handleId-1) / settings.plotGrid(2))+1, settings.plotGrid(1)) ~= 0);
          omitYLabel = (mod(handleId, settings.plotGrid(2)) ~= 1);
        else
          omitXLabel = false;
          omitYLabel = false;
        end
        % plot results
        isNotEmptyData = onePlot(relativeData, f, d, settings, ...
                               actualDisp, handleId*splitLegend, omitYLabel, omitXLabel);
        % check if any and which data were plotted in at least one function
        % and dimension
        plottedInAny = plottedInAny | isNotEmptyData;
        % increase handleId if empty plots are omitted
        if settings.omitEmptyPlots && any(isNotEmptyData)
          handleId = handleId + 1;
        % otherwise delete existing figure
        elseif settings.omitEmptyPlots && ~any(isNotEmptyData)
          close(handle{handleId})
        end
      end
    end
    if any(strcmp(settings.legendOption, {'out', 'manyout'})) && any(plottedInAny)
      % how many data are plotted
      nToPlot = sum(plottedInAny);
      % maximal number of data in one legend
      if strcmp(settings.legendOption, 'out')
        maxNamesLegend = nToPlot;
      else
        maxNamesLegend = 16;
      end
      % divide names to necessery sets
      nLegends = ceil(nToPlot/maxNamesLegend);
      setNumbers = floor(nToPlot/nLegends)*ones(1, nLegends);
      remNumber = mod(nToPlot, nLegends);
      setNumbers(1:remNumber) = setNumbers(1:remNumber) + 1;
      % create legend id vector
      setBounds = [0, cumsum(setNumbers)];
      idToPlot = inverseIndex(plottedInAny);
      for l = 1:nLegends
        actualID = idToPlot(setBounds(l) + 1 : setBounds(l+1));
        handle{nPlots + l} = soloLegend(settings.colors(actualID, :), settings.datanames(actualID), 2);
      end
    end
  end
  
end

function notEmptyData = onePlot(relativeData, fId, dId, ...
                 settings, dispLegend, splitLegendStatus, omitYLabel, omitXLabel)
% Plots one scaled graph 
%
% Note: Omitting y-label is currently enabled. To change this status
% uncomment rows at the end of onePlot function.

  nRelativeData = length(relativeData);
  if (~exist('omitXLabel', 'var') || isempty(omitXLabel))
    omitXLabel = false;
  end

  % parsing settings
  maxEval = settings.maxEval;
  colors = settings.colors;
  datanames = settings.datanames;
  aggFuns = settings.aggFuns;
  aggDims = settings.aggDims;
  funcNames = settings.funcNames;
  BBfunc = settings.BBfunc;
  dims = settings.dims;
  lineSpec = settings.lineSpec;
  medianLineWidth = settings.lineWidth;
  if length(medianLineWidth) < nRelativeData
    medianLineWidth = medianLineWidth(1)*ones(1, nRelativeData);
  end
  if (any(settings.drawQuantiles))
    % we should draw 1st and 3rd quantiles, too, with the same
    % colors and thinner lines
    medianLineWidth(((end/3)+1):end) = ceil(0.5*medianLineWidth(((end/3)+1):end));
    % % This has been already solved:
    % colors = repmat(colors, 3, 1);
    % lineSpec = repmat(lineSpec, 1, 3);
  end
  fullLegend = any(strcmp(settings.legendOption, {'first', 'split'}));

  notEmptyData = false(1, nRelativeData);
  % find not empty data
  for dat = 1:nRelativeData
    notEmptyData(dat) = ~isempty(relativeData{dat}) && size(relativeData{dat}, 1) >= fId ...
        && size(relativeData{dat}, 2) >= dId && ~isempty(relativeData{dat}{fId, dId});
  end
  % if (any(settings.drawQuantiles))
  %   nDataOriginal = length(settings.drawQuantiles);
  %   notEmptyData = notEmptyData & [true(1, nDataOriginal), settings.drawQuantiles, settings.drawQuantiles];
  % end
  
  if any(notEmptyData)
    nEmptyId = inverseIndex(notEmptyData);
    nUsefulData = sum(notEmptyData);
    % find minimal number of function evaluations to plot
    evaldim = 1:min([arrayfun(@(x) length(relativeData{nEmptyId(x)}{fId, dId}), 1:nUsefulData), maxEval]);
    h = zeros(1, nRelativeData);
    N_MARKERS = 3;
    % plot the lines
    thisYMin = Inf;
    thisYMax = -Inf;
    for dat = nRelativeData:-1:1
      if (notEmptyData(dat) && (~any(settings.drawQuantiles) ...
          || (dat <= (nRelativeData/3)) || settings.drawQuantiles(mod(dat-1,ceil(nRelativeData/3))+1)))
        h(dat) = plot(evaldim, relativeData{dat}{fId, dId}(evaldim), ...
          lineSpec{dat}, 'LineWidth', medianLineWidth(dat), 'Color', colors(dat, :));
        thisY = relativeData{dat}{fId, dId}(evaldim(:));
        thisY = thisY(:);
        thisYMin = min([thisYMin; thisY]);
        thisYMax = max([thisYMax; thisY]);
      else
        h(dat) = plot(0, 0, ...
          lineSpec{dat}, 'LineWidth', medianLineWidth(dat), 'Color', colors(dat, :), ...
          'Visible', 'off');
      end
      if (dat == nRelativeData)
        hold on;
        grid on;
      end
      if (~isempty(settings.markers) && (length(settings.markers) >= dat) ...
          && (ischar(settings.markers{dat}) && ~strcmpi(settings.markers{dat}, '')))
        randomStartMarker = 1+floor(rand()*(evaldim(end) / (N_MARKERS+1)));
        randomStepMarker  = floor((0.7 + rand()*0.6) * (evaldim(end) / (N_MARKERS+1)));
        randomXMarkers    = randomStartMarker:randomStepMarker:evaldim(end);
        if notEmptyData(dat)
          h(dat) = plot(randomXMarkers, relativeData{dat}{fId, dId}(randomXMarkers), ...
            settings.markers{dat}, 'LineWidth', medianLineWidth(dat), 'Color', colors(dat, :));
        else
          h(dat) = plot(0, 0, ...
            settings.markers{dat}, 'LineWidth', medianLineWidth(dat), 'Color', colors(dat, :), ...
            'Visible', 'off');
        end
      end
    end
    % if the legend should be plotted with all labels given, add artificial
    % not visible data
    if fullLegend && any(~notEmptyData)
      legendData = true(1, nRelativeData);
    else
      legendData = notEmptyData;
    end
    if (any(settings.drawQuantiles))
      % nRelativeData = nRelativeData / 3;
      legendData(((end/3)+1):end) = false;
    end
    nLegendData = sum(legendData);
    % legend settings
    if dispLegend
      legIds = false(size(legendData));
      switch splitLegendStatus
        % full legend
        case 0
          legIds = legendData;
        % first half of the legend
        case 1
          legIds(find(legendData, ceil(nLegendData/2), 'first')) = true;
        % second half of the legend
        case 2
          legIds(find(legendData, nLegendData - ceil(nLegendData/2), 'last')) = true;
      end
      if any(legIds)
        legend(h(legIds), datanames(legIds), 'Location', settings.legendLocation)
      end
    end
    ax = gca();
    if (settings.scaleY08)
      ax.YLim = [-8, 0];
    else
      ax.YLim = [thisYMin, thisYMax];
    end
  else
    warning('Function %d dimension %d has no data available.', BBfunc(fId), dims(dId))
  end

  % make title
  titleString = defopts(settings, 'TitleString', '');
  if isempty(titleString)
      if ~aggFuns
          titleString = ['f', num2str(BBfunc(fId))];
          if funcNames
              titleString = [titleString, ' ', bbobFunctionName(BBfunc(fId))];
          end
      end
      if ~aggDims
          titleString = [titleString, ' ', num2str(dims(dId)),'-D'];
      end
  end
  title(titleString)
  if (~omitXLabel)
    xlabel('Počet ohodnocení/ n');
  end
  if (~omitYLabel)
    ylabel('\Delta_f^{log}');
  end
  
  hold off
  
  if (any(settings.drawQuantiles))
    notEmptyData = notEmptyData(1:(end/3));
  end
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
  switch nNames
    case 1
      fb = [3.5, -2]; % [fb_x, fb_y]
    case 2
      fb = [1.75, -1.5]; % [fb_x, fb_y]
    otherwise
      fb = [1, -2*(nNames+1)*yDiff]; % [fb_x, fb_y]
  end
  line([fa(1), fb(1)], [fa(2), fa(2)], 'Color', 'k') %  _
  line([fa(1), fa(1)], [fa(2), fb(2)], 'Color', 'k') % |
  line([fb(1), fb(1)], [fa(2), fb(2)], 'Color', 'k') %  _|
  line([fa(1), fb(1)], [fb(2), fb(2)], 'Color', 'k') %  
  
  hold off
end

function actualDisp = dispLegend(plotId, legendOption)
% returns if the legend should be displayed

  % do not display legend inside the graph in 'hide', 'out', 'manyout' 
  % cases
  if any(strcmp(legendOption, {'hide', 'out', 'manyout'}))
    displayLegend = false;
  else
    displayLegend = true;
  end
  
  % figure 1  - 'show', 'split', 'first'
  % figure 2  - 'show', 'split'
  % figure 3+ - 'show'
  actualDisp = displayLegend && ...
               ( strcmp(legendOption, 'show') || ...
                (plotId == 1) || ...
                (plotId == 2 && strcmp(legendOption, 'split')) ...
               );  
end

function name = bbobFunctionName(fId)
% returns name of function fId
  switch fId
    % separable functions
    case 1
      name = 'Sphere';
    case 2
      name = 'Ellipsoidal';
    case 3 	
      name = 'Rastrigin';
    case 4
      name = 'Bueche-Rastrigin';
    case 5
      name = 'Linear Slope';
      
    % functions with low or moderate conditioning
    case 6
      name = 'Attractive Sector';
    case 7
      name = 'Step Ellipsoidal';
    case 8
      name = 'Rosenbrock, original';
    case 9 
      name = 'Rosenbrock, rotated';
      
    % functions with high conditioning and unimodal
    case 10
      name = 'Ellipsoidal, high conditioning';
    case 11
      name = 'Discus';
    case 12
      name = 'Bent Cigar';
    case 13
      name = 'Sharp Ridge';
    case 14
      name = 'Different Powers';
      
    % multi-modal functions with adequate global structure
    case 15
      name = 'Rastrigin, multi-modal';
    case 16
      name = 'Weierstrass';
    case 17
      name = 'Schaffers F7';
    case 18
      name = 'Schaffers F7, ill-conditioned';
    case 19
      name = 'Composite Griewank-Rosenbrock F8F2';
    
    % multi-modal functions with weak global structure
    case 20
      name = 'Schwefel';
    case 21
      name = 'Gallagher''s Gaussian 101-me Peaks';
    case 22
      name = 'Gallagher''s Gaussian 21-hi Peaks';
    case 23
      name = 'Katsuura';
    case 24
      name = 'Lunacek bi-Rastrigin';
      
    % TODO: noisy functions
    % other functions
    otherwise
      name = '';
  end
end
