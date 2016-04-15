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
  minValue = defopts(settings, 'MinValue', 10^(-8));
  maxEval = defopts(settings, 'MaxEval', 250);
  statistic = defopts(settings, 'Statistic', @mean);
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
  data_stats = cellfun(@(D) gainStatistic(D, dimIds, funcIds, useMaxInstances, false, statistic), ...
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
  handle = relativePlot(data_stats, dims, BBfunc, length(funcIds), datanames, colors, aggDims, aggFuns, maxEval);

end

function handle = relativePlot(data_stats, dims, BBfunc, numOfFuncIds, datanames, colors, aggDims, aggFuns, maxEval)
% Plots quantile graph of different algorithms in one function and one
% dimension

  numOfData = length(data_stats);
  evaldim = 1:length(data_stats{1}{1});
  medianLineWidth = 2;
  minGraph = 10e-8;
  maxGraph = 1;
  
  for f = 1:numOfFuncIds
    % find useful data and plot 
    for d = 1:length(dims)
      % find available data
      notEmptyData = true(1, numOfData);
      for dat = 1:numOfData
        notEmptyData(dat) = ~isempty(data_stats{dat}{f,d});
      end
      if any(notEmptyData)
        nEmptyId = inverseIndex(notEmptyData);
        nUsefulData = sum(notEmptyData);
        % count f-values ratio
        actualData = cell2mat(arrayfun(@(D) data_stats{nEmptyId(D)}{f,d}, 1:nUsefulData, 'UniformOutput', false));
        actualMin = min(min(actualData));
        actualMax = max(max(actualData));
        for D = 1:nUsefulData
          relativeData{nEmptyId(D)}{f, d} = log(((actualData(:, D) - actualMin) * (maxGraph - minGraph)/(actualMax - actualMin))' + minGraph);
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
        relativeData{D}{1, d} = mean(cell2mat(arrayfun(@(x) relativeData{D}{x,d}, 1:nFuns, 'UniformOutput', false)'));
      end
      relativeData{D} = relativeData{D}(1, :);
    end
  else
    nFunsToPlot = length(BBfunc);
  end
  
  handle = zeros(1, nDimsToPlot*nFunsToPlot);
  for f = 1:nFunsToPlot
    for d = 1:nDimsToPlot
      handle((d-1) * nFunsToPlot + f) = figure('Units', 'centimeters', 'Position', [1 1 12.5 6]);
      for dat = 1:numOfData
        notEmptyData(dat) = ~isempty(relativeData{dat}{f, d});
      end
      if any(notEmptyData)
        h = zeros(1, nUsefulData);
        ftitle = cell(1, nUsefulData);
        % mean
        h(1) = plot(evaldim, relativeData{nEmptyId(1)}{f, d}(1:maxEval), ...
          'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(1), :));
        ftitle{1} = datanames{nEmptyId(1)};
        hold on
        grid on
        for dat = 2:nUsefulData
          h(dat) = plot(evaldim, relativeData{nEmptyId(dat)}{f, d}(1:maxEval), ...
            'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(dat), :));
          ftitle{dat} = datanames{nEmptyId(dat)};
        end
        legend(h, ftitle, 'Location', 'NorthEast')
      else
        warning('Function %d dimension %d has no data available', BBfunc(f), dims(d))
      end
      
      titleString = '';
      if ~aggFuns
        titleString = ['f', num2str(BBfunc(f))];
      end
      if ~aggDims
        titleString = [titleString, ' ', num2str(dims(d)),'D'];
      end
      title(titleString)
      xlabel('Number of evaluations / D')
      ylabel('Minimum function values')
      hold off
    end
  end
end