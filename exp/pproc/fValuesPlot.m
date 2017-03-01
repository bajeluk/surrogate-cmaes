function handle = fValuesPlot(data, varargin)
% handle = fValuesPlot(data, settings)
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
%     'Dependency'    - 'algorithm' or 'dimension'
%
% Output:
%   handle - handles of resulting figures
%
% See Also:
%   speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  if nargin < 1 || isempty(data)
    handle = [];
    help fValuesPlot
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
  avgDims = defopts(settings, 'AggregateDims', false);
  dependency = defopts(settings, 'Dependency', 'algorithms');
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
  data_stats = cellfun(@(D) gainStatistic(D, dimIds, funcIds, ...
                                          'MaxInstances', useMaxInstances, ...
                                          'AverageDims', avgDims, ...
                                          'Statistic', statistic), ...
                            data, 'UniformOutput', false);

  switch lower(dependency)
    case {'alg', 'algorithm'}
      handle = quantileAlgorithmPlot(data_stats, dims, BBfunc, length(funcIds), datanames, colors, avgDims);
    case {'dim', 'dimension'}
      handle = quantileDimensionPlot(data_stats, dims, BBfunc, length(funcIds), datanames, colors);
  end

end

function handle = quantileAlgorithmPlot(data_stats, dims, BBfunc, numOfFuncIds, datanames, colors, avgDims)
% Plots quantile graph of different algorithms in one function and one
% dimension

  numOfData = length(data_stats);
  evaldim = 1:length(data_stats{1}{1});
  medianLineWidth = 2;
  doQuantiles = size(data_stats{1}{1,1}, 2) == 3;
  if doQuantiles
    quantileLineWidth = 0.5;
    quantileStyle = '-';
    mainLine = 2;
  else
    mainLine = 1;
  end
  
  if avgDims
    nDimsToPlot = 1;
  else
    nDimsToPlot = length(dims);
  end
  
  handle = zeros(1, nDimsToPlot*numOfFuncIds);
  for d = 1:nDimsToPlot
    for f = 1:numOfFuncIds
      handle((d-1) * numOfFuncIds + f) = figure('Units', 'centimeters', 'Position', [1 1 12.5 6]);
      % find available data
      notEmptyData = true(1, numOfData);
      for dat = 1:numOfData
        notEmptyData(dat) = ~isempty(data_stats{dat}{f,d});
      end
      if any(notEmptyData)
        nEmptyId = inverseIndex(notEmptyData);
        h = zeros(1, sum(notEmptyData));
        ftitle = cell(1, sum(notEmptyData));
        % median
        h(1) = semilogy(evaldim, data_stats{nEmptyId(1)}{f, d}(evaldim, mainLine), ...
          'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(1), :));
        ftitle{1} = datanames{nEmptyId(1)};
        hold on
        grid on
        % quantiles
        if doQuantiles
          semilogy(evaldim, data_stats{nEmptyId(1)}{f, d}(evaldim, 1), ...
            'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(1), :));
          semilogy(evaldim, data_stats{nEmptyId(1)}{f, d}(evaldim, 3), ...
            'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(1), :));
        end
        for dat = 2:sum(notEmptyData)
          h(dat) = semilogy(evaldim, data_stats{nEmptyId(dat)}{f, d}(evaldim, mainLine), ...
            'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(dat), :));
          ftitle{dat} = datanames{nEmptyId(dat)};
          if doQuantiles
            semilogy(evaldim, data_stats{nEmptyId(dat)}{f, d}(evaldim, 1), ...
              'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(dat), :));
            semilogy(evaldim, data_stats{nEmptyId(dat)}{f, d}(evaldim, 3), ...
              'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(dat), :));
          end
        end

        % additional plot settings
        % ylim(gca,[1e-8 1e5])

        legend(h, ftitle, 'Location', 'NorthEast')
      else
        warning('Function %d dimension %d has no data available', BBfunc(f), dims(d))
      end
      
      if avgDims
        title(['f', num2str(BBfunc(f))])
      else
        title(['f', num2str(BBfunc(f)), ' ', num2str(dims(d)),'D'])
      end
      xlabel('Number of evaluations / D')
      ylabel('Minimum function values')
      hold off
    end
  end
end

function handle = quantileDimensionPlot(data_stats, dims, BBfunc, numOfFuncIds, datanames, colors)
% Plots quantile graph of different dimensions in one function and one data

  numOfData = length(data_stats);
  numOfDims = length(dims);
  evaldim = 1:length(data_stats{1}{1});
  medianLineWidth = 2;
  doQuantiles = size(data_stats{1}{1,1}, 2) == 3;
  if doQuantiles
    quantileLineWidth = 0.5;
    quantileStyle = '-';
    mainLine = 2;
  else
    mainLine = 1;
  end
  
  handle = zeros(1, numOfData*numOfFuncIds);
  for dat = 1:numOfData
    for f = 1:numOfFuncIds
      handle((dat-1) * numOfFuncIds + f) = figure('Units', 'centimeters', 'Position', [1 1 12.5 6]);
      % find available data
      notEmptyDim = true(1, numOfDims);
      for d = 1:numOfDims
        notEmptyDim(d) = ~isempty(data_stats{dat}{f,d});
      end
      if any(notEmptyDim)
        nEmptyId = inverseIndex(notEmptyDim);
        h = zeros(1, sum(notEmptyDim));
        ftitle = cell(1, sum(notEmptyDim));
        % median
        h(1) = semilogy(evaldim, data_stats{dat}{f, nEmptyId(1)}(evaldim, mainLine), ...
          'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(1), :));
        ftitle{1} = [num2str(dims(nEmptyId(1))), 'D'];
        hold on
        grid on
        % quantiles
        if doQuantiles
          semilogy(evaldim, data_stats{dat}{f, nEmptyId(1)}(evaldim, 1), ...
            'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(1), :));
          semilogy(evaldim, data_stats{dat}{f, nEmptyId(1)}(evaldim, 3), ...
            'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(1), :));
        end
        for d = 2:sum(notEmptyDim)
          h(d) = semilogy(evaldim, data_stats{dat}{f, nEmptyId(d)}(evaldim, mainLine), ...
            'LineWidth', medianLineWidth, 'Color', colors(nEmptyId(d), :));
          ftitle{d} = [num2str(dims(nEmptyId(d))), 'D'];
          if doQuantiles
            semilogy(evaldim, data_stats{dat}{f, nEmptyId(d)}(evaldim, 1), ...
              'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(d), :));
            semilogy(evaldim, data_stats{dat}{f, nEmptyId(d)}(evaldim, 3), ...
              'LineStyle', quantileStyle, 'LineWidth', quantileLineWidth, 'Color', colors(nEmptyId(d), :));
          end
        end

        % additional plot settings
        % ylim(gca,[1e-8 1e5])

        legend(h, ftitle, 'Location', 'NorthEast')
      else
        warning('Function %d algorithm %s has no data available', BBfunc(f), datanames{dat})
      end
      
      title([datanames{dat}, ' f', num2str(BBfunc(f))])
      xlabel('Number of evaluations / D')
      ylabel('Minimum function values')
      hold off
    end
  end
end