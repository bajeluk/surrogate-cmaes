function handle = metafeaturePlot(mftsVal, dataVal, varargin)
% handle = metafeaturePlot(mftsVal, errVal, settings)
% Plot dependency of model error on metafeature values.
%
% Input:
%   data     - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataColor'     - vector Nx3 of (RGB) colors of individual 
%                       algorithms, N is the number of algorithms | double 
%                       array
%                       (e.g. 2 algorithms: [115, 60, 215; 0, 40, 60]) 
%     'DataNames'     - cell array of data names (e.g. names of algorithms)
%                       | cell array of strings
%                       (e.g. {'alg1', 'alg2', 'alg3'})
%     'MedianLS'      - median line style (specification) | double
%     'MedianLW'      - median line width | double
%     'NValues'       - number of values to plot | integer
%     'QuartileLS'    - quartile line style (specification) | double
%     'QuartileLW'    - quartile line width | double
%     'ShowLegend'    - show plot legend | boolean
%
% Output:
%   handle - handles of resulting figures | cell-array
%
% See Also:
%   relativeFValuesPlot


  % initialization
  if nargin < 2 || isempty(mftsVal) || isempty(dataVal)
    if nargout > 0
      handle = {};
    end
    help relativeFValuesPlot
    return
  end
  settings = settings2struct(varargin);

  % parse metafeatures
  nFeat = size(mftsVal, 2);
  if istable(mftsVal)
    defMftsNames = mftsVal.Properties.VariableNames;
    mftsVal = table2array(mftsVal);
  else
    defMftsNames = arrayfun(@(x) ['mfts_', num2str(x)], 1:nFeat, 'UniformOutput', false);
  end
  % parse error values
  nData = size(dataVal, 2);
  if istable(dataVal)
    defDataNames = dataVal.Properties.VariableNames;
    dataVal = table2array(dataVal);
  else
    defDataNames = arrayfun(@(x) ['data_', num2str(x)], 1:nData, 'UniformOutput', false);
  end
  
  % parse settings
  kerColor = defopts(settings, 'DataColor', getAlgColors(1:nData)/255);
  dataNames = defopts(settings, 'DataNames', defDataNames);
  medianLineWidth = defopts(settings, 'MedianLW', 1.8);
  quartileLineWidth = defopts(settings, 'QuartileLW', 1);
  medianLineStyle = defopts(settings, 'MedianLS', '-');
  quartileLineStyle = defopts(settings, 'QuartileLS', '-.');
  nPointsToPlot = defopts(settings, 'NValues', 200);
  logBound = defopts(settings, 'LogBound', 5);
  mftsNames = defopts(settings, 'MftsNames', defMftsNames);
  showLegend = defopts(settings, 'ShowLegend', true);
  
  % check names format
  if ~iscell(dataNames)
    dataNames = {dataNames};
  end
  if ~iscell(mftsNames)
    mftsNames = {mftsNames};
  end
  assert(numel(dataNames) == nData, 'Number of DataNames is not equal to the number of data')
  assert(numel(mftsNames) == nFeat, 'Number of MftsNames is not equal to the number of features')
  
  handle = cell(1, nFeat);
  % feature loop
  for mf = 1:nFeat
    % get default nPointsToPlot value
    nPointsToPlotAct = nPointsToPlot;
    
    % remove metafeature NaN and Inf values
    nanOrInf = isnan(mftsVal(:, mf)) | isinf(mftsVal(:, mf));
    actual_mftsVal = mftsVal(~nanOrInf, mf);
    % create x-values for plot
    if numel(unique(actual_mftsVal)) < nPointsToPlotAct
      nPointsToPlotAct = numel(unique(actual_mftsVal));
      xBound = [min(actual_mftsVal) - 1, unique(actual_mftsVal')]; % quantile(actual_dataVal, nPointsToPlotAct)];
      logAxis = false;
    elseif abs(diff(log10(minmax(actual_mftsVal(~isinf(actual_mftsVal))')))) > logBound      
      % xBound = logspace(log10(min(actual_dataVal) - eps), log10(max(actual_dataVal)), nPointsToPlotAct + 1);
      xBound = [min(actual_mftsVal) - 1, quantile(actual_mftsVal, nPointsToPlotAct)];
      logAxis = true;
    else
      % xBound = linspace(min(actual_dataVal) - eps, max(actual_dataVal), nPointsToPlotAct + 1);
      xBound = [min(actual_mftsVal) - 1, quantile(actual_mftsVal, nPointsToPlotAct)];
      logAxis = false;
    end
    yVal = NaN(nPointsToPlotAct, 3);
  
    % name of first model error
    actual_dataVal = dataVal(~nanOrInf, 1);
    % create median and quartiles of data values for plot
    for xb = 1:nPointsToPlotAct
      yVal(xb, :) = quantile((actual_dataVal(actual_mftsVal >  xBound(xb) & ...
                                             actual_mftsVal <= xBound(xb+1))), ...
                             [1/4, 1/2, 3/4]);
    end
    % check missing values
    plotPointId = ~all(isnan(yVal), 2);
    
    % create figure
    handle{mf} = figure();
    h = [];
    if logAxis
      h(1) = semilogx(xBound([false; plotPointId]), yVal(plotPointId, 1), ...
        'Color', kerColor(1, :), ...
        'LineStyle', quartileLineStyle, ...
        'LineWidth', quartileLineWidth);
      hold on
      h(end+1) = semilogx(xBound([false; plotPointId]), yVal(plotPointId, 2), ...
        'Color', kerColor(1, :), ...
        'LineStyle', medianLineStyle, ...
        'LineWidth', medianLineWidth);
      h(end+1) = semilogx(xBound([false; plotPointId]), yVal(plotPointId, 3), ...
        'Color', kerColor(1, :), ...
        'LineStyle', quartileLineStyle, ...
        'LineWidth', quartileLineWidth);
    else
      h(1) = plot(xBound([false; plotPointId]), yVal(plotPointId, 1), ...
        'Color', kerColor(1, :), ...
        'LineStyle', quartileLineStyle, ...
        'LineWidth', quartileLineWidth);
      hold on
      h(end+1) = plot(xBound([false; plotPointId]), yVal(plotPointId, 2), ...
        'Color', kerColor(1, :), ...
        'LineStyle', medianLineStyle, ...
        'LineWidth', medianLineWidth);
      h(end+1) = plot(xBound([false; plotPointId]), yVal(plotPointId, 3), ...
        'Color', kerColor(1, :), ...
        'LineStyle', quartileLineStyle, ...
        'LineWidth', quartileLineWidth);
    end
  
    % add the rest of models
    for m = 2:nData
      actual_dataVal = dataVal(~nanOrInf, m);
  
      % create median and quartiles of data values for plot
      for xb = 1:nPointsToPlotAct
        yVal(xb, :) = quantile((actual_dataVal(actual_mftsVal > xBound(xb) & actual_mftsVal <= xBound(xb+1))), [1/4, 1/2, 3/4]);
      end
      % check missing values
      plotPointId = ~all(isnan(yVal), 2);
      if logAxis
        h(end+1) = semilogx(xBound([false; plotPointId]), yVal(plotPointId, 1), ...
          'Color', kerColor(m, :), ...
          'LineStyle', quartileLineStyle, ...
          'LineWidth', quartileLineWidth);
        h(end+1) = semilogx(xBound([false; plotPointId]), yVal(plotPointId, 2), ...
          'Color', kerColor(m, :), ...
          'LineStyle', medianLineStyle, ...
          'LineWidth', medianLineWidth);
        h(end+1) = semilogx(xBound([false; plotPointId]), yVal(plotPointId, 3), ...
          'Color', kerColor(m, :), ...
          'LineStyle', quartileLineStyle, ...
          'LineWidth', quartileLineWidth);
      else
        h(end+1) = plot(xBound([false; plotPointId]), yVal(plotPointId, 1), ...
          'Color', kerColor(m, :), ...
          'LineStyle', quartileLineStyle, ...
          'LineWidth', quartileLineWidth);
        h(end+1) = plot(xBound([false; plotPointId]), yVal(plotPointId, 2), ...
          'Color', kerColor(m, :), ...
          'LineStyle', medianLineStyle, ...
          'LineWidth', medianLineWidth);
        h(end+1) = plot(xBound([false; plotPointId]), yVal(plotPointId, 3), ...
          'Color', kerColor(m, :), ...
          'LineStyle', quartileLineStyle, ...
          'LineWidth', quartileLineWidth);
      end
    end
    title(strrep(mftsNames{mf}, '_', '\_'))
    xlabel('feature value')
    ylabel('model RDE')
    if showLegend
      legend(h(2:3:end), dataNames, 'Location', 'Best')
    end
    axis([xBound(2), xBound(end), 0, 1])
    hold off
  
  end
  
end

function h = plotLine(plotType, xBound, yVal, plotPointId)
% plot median or quartile line
  h = plot(xBound([false; plotPointId]), yVal(plotPointId, 3), ...
           'Color', kerColor(m, :), ...
           'LineStyle', quartileLineStyle, ...
           'LineWidth', quartileLineWidth);
end