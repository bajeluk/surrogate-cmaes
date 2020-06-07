function mftPropPlot(value, nans, inf_plus, inf_mins, varargin)
% mftPropPlot - plot metafeature property
%
% mftPropPlot(value, nans, inf_plus, inf_mins, settings)
%
% Input:
%   value    - metafeature property values | double vector nObservations X
%              nDifferentValues to plot
%   nans     - number of nan values | integer vector
%   inf_plus - number of Inf values | integer vector
%   inf_mins - number of -Inf values | integer vector
%   settings - pairs of property (string) and value, or struct with 
%              properties as fields:
%     'MftName'   - metafeature name | string
%     'ValueName' - name of property value | string

  if nargin < 4
    if nargin < 3
      if nargin < 2
        if nargin < 1
          help mftPropPlot
          return
        end
        nans = [];
      end
      inf_plus = [];
    end
    inf_mins = [];
  end
  
  % init variables
  [nObs, nLines] = size(value);
  plotVec = 1:nObs;
  plotNan = any(nans > 0);
  plotInfP = any(inf_plus > 0);
  plotInfM = any(inf_mins > 0);
  
  % define colors
  valueColor = [0, 0.4470, 0.7410];
  nanColor = [0.8500, 0.3250, 0.0980];
  infColor = [0.4660, 0.6740, 0.1880];
  dimColor = [233, 188, 213]/256;
  dimLabelColor = [195, 72, 141]/256;
  
  % parse settings
  settings = settings2struct(varargin{:});
  dimensions = defoptsi(settings, 'Dimensions', []);
  dimTypes  = defoptsi(settings, 'DimTypes', []);
  mftName = defoptsi(settings, 'MftName', 'metafeature');
  defValueNames = compose('value%d', 1:nLines);
  valueNames = defoptsi(settings, 'ValueName', defValueNames);
  valueStyle = defoptsi(settings, 'ValueStyle', repmat({'-'}, 1, nLines));
  xValues = defoptsi(settings, 'XValues', plotVec);
  xlabelName = defoptsi(settings, 'XLabelName', 'quantile');
  
  nDims = numel(dimTypes);
  
  % check value names -> cell-array of nLines lenght
  if ~iscell(valueNames) && ischar(valueNames)
    valueNames = {valueNames};
  elseif ~iscell(valueNames)
    warning('ValueNames are not properly defined strings. Using defaults.')
    valueNames = defValueNames;
  end
  if numel(valueNames) < nLines
    warning('Not enough value names. Filling the rest with defaults.')
    valueNames(end+1:nLines) = defValueNames(numel(valueNames)+1 : end);
  end

  % display dimensions indicator
  plotDims = ~isempty(dimensions) && ( nDims == size(dimensions, 2) );
  
  % axis ids
  if plotDims
    valueAxisId = 2;
    nanInfAxisId = 1;
  else
    valueAxisId = 1;
    nanInfAxisId = 2;
  end
  
  % create plot
  figure()
  title(mftName)
  hax = gca;
  
  % line handle
  h = [];
  % legend list
  legList = {};
  legListNI = {};
  
  hold on  
  % add dimension ranges markers
  if plotDims
    yyaxis left
    
    plotDimensionsGrid = [];
    for d = 1:nDims
      plotDimensionsGrid = [plotDimensionsGrid; repmat(dimensions(:, d)', 20, 1)];
    end
    imagesc(plotDimensionsGrid)
    hold on
    colormap([1,1,1; dimColor])
    % return own axis
    line([0, 1000], [0, 0], 'Color', [0, 0, 0])       %       x-axis
    line([1000, 1000], [0, 100], 'Color', valueColor) % right y-axis
    if plotNan
      line([0, 0], [0, 100], 'Color', nanColor)       %  left y-axis
    end
    if plotInfP || plotInfM
      line([0, 0], [0, 100], 'Color', infColor)       %  left y-axis
    end
    % label dimensions
    dimLabel = strsplit(num2str(dimTypes, '%dD,'), ',');
    dimLabel = dimLabel(1:end-1); % last cell is empty
    if plotNan || plotInfP || plotInfM
      dimXLabelLocation = 50*ones(1, nDims);
    else
      dimXLabelLocation = -100*ones(1, nDims);
    end
    dimYLabelLocation = 100/nDims*(1:nDims) - 50/nDims;
    text(dimXLabelLocation, dimYLabelLocation, dimLabel, ...
         'Color', dimLabelColor, 'FontSize', 14)
    % normalize to [0, 100] range
    hax.YAxis(1).Limits = [0, 100];
  elseif ( nDims ~= size(dimensions, 2) )
    warning('Dimension columns and number of dimTypes differ. Ignoring dimension settings.')
  else
    
  end
  
  % plot values
  if plotDims
    yyaxis right
  else
    yyaxis left
  end
  hax.YAxis(valueAxisId).Color = valueColor;
  % plot individual lines
  for l = 1:nLines
    h(end+1) = plot(plotVec, value(plotVec, l), ...
                'Color', valueColor, ...
                'LineStyle', valueStyle{l});
    legList(end+1) = valueNames(l);
  end
  % check resulting ylimits
  hax.YAxis(valueAxisId).Limits = [max(hax.YAxis(valueAxisId).Limits(1), 0), ...
                                   min(hax.YAxis(valueAxisId).Limits(2), 1)];
  hax.XAxis.Limits = [0, max(plotVec)];
  
  % title and axis labels
  xlabel(xlabelName)
  ylabel('Feature value')
  
  % plot NaN and Inf values
  if plotNan || plotInfP || plotInfM
    % create extra axis
    if plotDims
      yyaxis left
    else
      yyaxis right
    end
    % normalize to [0, 100] range
    hax.YAxis(nanInfAxisId).Limits = [0, 100];
    yLabelList = {};
    legListNI = {};
    
    % plot NaNs
    if plotNan
      legListNI(end+1) = {'0.5 NaN'};
      h(end+1) = area(plotVec, nans, 'LineWidth', 1, ...
                                      'FaceAlpha', 0.2, ...
                                      'FaceColor', nanColor, ...
                                      'EdgeColor', nanColor ...
                     );
      yLabelList{end+1} = 'NaN';
      hax.YAxis(nanInfAxisId).Color = nanColor;
    end
    % plot Infs
    if plotInfP
      legListNI(end+1) = {'0.5 Inf'};
      h(end+1) = area(plotVec, inf_plus, 'LineWidth', 1, ...
                                         'FaceColor', infColor, ...
                                         'EdgeColor', infColor, ...
                                         'FaceAlpha', 0.2);
    end
    % plot -Infs
    if plotInfM
      legListNI(end+1) = {'0.5 -Inf'};
      h(end+1) = area(plotVec, inf_mins, 'LineWidth', 1, ...
                                         'FaceColor', infColor, ...
                                         'EdgeColor', infColor, ...
                                         'FaceAlpha', 0.2);
    end
    if plotInfP || plotInfM
      yLabelList{end+1} = 'Inf';
      hax.YAxis(nanInfAxisId).Color = infColor;
    end
    
    % ylabel
    ylabel(sprintf('%% of %s samples', strjoin(yLabelList, ' or ')))
  else
    % disable axis for NaN or Inf
    hax.YAxis(nanInfAxisId).Visible = 'off';
  end
  
  % show legend
  legend(h, [legList, legListNI])
  
  hold off

end