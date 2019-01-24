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
%     'LogBound'      - exponential bound for logarithmic scale | double
%     'MedianLS'      - median line style (specification) | double
%     'MedianLW'      - median line width | double
%     'MftsNames'     - metafeature names | cell-array of string
%     'NValues'       - number of values to plot | integer
%     'QuantileRange' - range of metafeatures according to quantile values 
%                       | double (e.g. 5% interval [0.05, 0.95])
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
  plotSet.medianLineWidth = defopts(settings, 'MedianLW', 1.8);
  plotSet.quartileLineWidth = defopts(settings, 'QuartileLW', 1);
  plotSet.medianLineStyle = defopts(settings, 'MedianLS', '-');
  plotSet.quartileLineStyle = defopts(settings, 'QuartileLS', '-.');
  nPointsToPlot = defopts(settings, 'NValues', 200);
  logBound = defopts(settings, 'LogBound', 5);
  mftsNames = defopts(settings, 'MftsNames', defMftsNames);
  quantRange = defopts(settings, 'QuantileRange', [0, 1]);
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
    
    % remove metafeature values out of quantile range (automatically
    % removes Inf and NaN)
    quantBounds = quantile(mftsVal(:, mf), quantRange);
    inBounds = (mftsVal(:, mf) >= quantBounds(1) & ...
                mftsVal(:, mf) <= quantBounds(2));
    actual_mftsVal = mftsVal(inBounds, mf);
    % prepare minimal and maximal value
    mnmx = minmax(actual_mftsVal');
    % create x-values for plot
    % small number of values
    if numel(unique(actual_mftsVal)) < nPointsToPlotAct
      nPointsToPlotAct = numel(unique(actual_mftsVal));
      plotSet.xBound = [mnmx(1) - 1, unique(actual_mftsVal')]; % quantile(actual_dataVal, nPointsToPlotAct)];
      plotSet.logAxis = false;
    % logarithmic scale in case large exponential differences and no
    % positive and negative values at the same moment
    elseif abs(diff(log10(mnmx))) > logBound && ...
           abs(sum(sign(mnmx))) > 0
      % xBound = logspace(log10(min(actual_dataVal) - eps), log10(max(actual_dataVal)), nPointsToPlotAct + 1);
      plotSet.xBound = [mnmx(1) - 1, quantile(actual_mftsVal, nPointsToPlotAct)];
      plotSet.logAxis = true;
    % linear scale
    else
      % xBound = linspace(min(actual_dataVal) - eps, max(actual_dataVal), nPointsToPlotAct + 1);
      plotSet.xBound = [mnmx(1) - 1, quantile(actual_mftsVal, nPointsToPlotAct)];
      plotSet.logAxis = false;
    end
    plotSet.yVal = NaN(nPointsToPlotAct, 3);
    
    % create figure
    handle{mf} = figure();
    h = [];
  
    % data loop
    for m = 1:nData
      actual_dataVal = dataVal(inBounds, m);
  
      % create median and quartiles of data values for plot
      for xb = 1:nPointsToPlotAct
        plotSet.yVal(xb, :) = quantile(...
          (actual_dataVal(actual_mftsVal >  plotSet.xBound(xb) & ...
                          actual_mftsVal <= plotSet.xBound(xb+1))), ...
           [1/4, 1/2, 3/4]);
      end
      % check missing values
      plotSet.plotPointId = ~all(isnan(plotSet.yVal), 2);
      % prepare color
      plotSet.kerColor = kerColor(m, :);
      % plot data
      h(end+1) = plotLine(1, plotSet);
      hold on
      h(end+1) = plotLine(2, plotSet);
      h(end+1) = plotLine(3, plotSet);
    end
    
    % additional figure captions and ranges
    title(strrep(mftsNames{mf}, '_', '\_'))
    xlabel('feature value')
    ylabel('model RDE')
    if showLegend
      legend(h(2:3:end), dataNames, 'Location', 'Best')
    end
    axis([plotSet.xBound(2), plotSet.xBound(end), 0, 1])
    hold off
  
  end % mfts loop
  
end

function h = plotLine(quartNumber, plotSet)
% plot median or quartile line
  plotInput = {...
           plotSet.xBound([false; plotSet.plotPointId]), ...
           plotSet.yVal(plotSet.plotPointId, quartNumber), ...
           'Color', plotSet.kerColor ...
           };
  % median
  if quartNumber == 2
    plotInput = [plotInput, {...
                              'LineStyle', plotSet.medianLineStyle, ...
                              'LineWidth', plotSet.medianLineWidth ...
                            }];
  % quartiles
  else
    plotInput = [plotInput, {...
                              'LineStyle', plotSet.quartileLineStyle, ...
                              'LineWidth', plotSet.quartileLineWidth ...
                            }];
  end
  % logarithmic scale
  if plotSet.logAxis
    h = semilogx(plotInput{:});
  % linear scale
  else
    h = plot(plotInput{:});
  end
  
end