function h = footprintPlot(baseData, valueData, varargin)
% Plot data footprint:
% 
%   h = footprintPlot(baseData, valueData, settings) 
%
% Plot only showData footprint:
%
%   h = footprintPlot(baseData, showData, valueData, settings)
%
% Input:
%   baseData - all data coordinates in 2 component space | Nx2 double
%   showData - date coordinates in 2 component space colored according to 
%              valueData | Mx2 double
%             - if not set, showData = baseData (M = N)
%   valueData - values of showData | Mx1 double
%   settings  - pairs of property (string) and value, or struct with 
%               properties as fields:
%     'Title'     - resulting plot title | string
%     'Statistic' - statistic used to process values in points with 
%                   identical coordinates | function handle (default =
%                   @median)
%     'ValueName' - name of valueData | string 
%
% Output:
%   h - resulting figure | handle
% 
% See Also:
%   metaLearn_loadData


  if nargout > 0
    h = {};
  end
  if nargin < 1
    help footprintPlot
    return
  end
  
  if isvector(valueData)
    showData = baseData;
  else
    showData = valueData;
    valueData = varargin{1};
    varargin = varargin(2:end);
  end
  % parse settings
  settings = settings2struct(varargin{:});
  plotTitle = defoptsi(settings, 'title', '');
  barTitle = defoptsi(settings, 'valuename', '');
  stat = defoptsi(settings, 'statistic', @median);
  
  % find unique data and merge their values
  [showData, ~, uniId] = unique(showData, 'rows');
  mergedValues = NaN(size(showData, 1), 1);
  for i = 1:size(showData, 1)
    mergedValues(i) = stat(valueData(uniId == i));
  end
  
  % TODO: correct reduction according to percentile
  percOutId = mergedValues == min(mergedValues);
  mergedValues(percOutId) = [];
  showData(percOutId, :) = [];

  % init graphic values
  baseCol = [220, 220, 220]/255; % grey
  ptSize = 20;

  h = figure();
  % base data
  scatter(baseData(:, 1), baseData(:, 2), ptSize, baseCol, 'filled')
  hold on
  scatter(showData(:, 1), showData(:, 2), ptSize, mergedValues, 'filled')
  colormap jet
  hcb = colorbar;
  if ~isempty(barTitle)
    colorTitleHandle = get(hcb, 'Title');
    set(colorTitleHandle, 'String', barTitle);
  end
  title(plotTitle)
  xlabel('First PCA component')
  ylabel('Second PCA component')
  hold off
end