function handle = groupFValuesPlot(data, varargin)
% handle = groupFValuesPlot(data, settings)
% Plots dependences of minimal function values on function 
% evaluations / dimension for individual functions.
%
% Input:
%   data     - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'AggregateDims' - aggregate dimensions in plots | boolean
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
    help groupFValuesPlot
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
  [nInputFuns, nInputDims] = size(data{1});
  numOfData = numel(data);
  defaultDims = [2, 3, 5*2.^(0:nInputDims-3)];
  dims = defopts(settings, 'DataDims', defaultDims(1:nInputDims));
  defaultFuns = [1:24, 101:130, 200 + (1:nInputFuns-54)];
  funs = defopts(settings, 'PlotFuns', defaultFuns(1:nInputFuns));
  % name settings
  datanames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:numOfData, 'UniformOutput', false));
  if length(datanames) ~= numOfData && any(emptyData)
    datanames = datanames(~emptyData);
  end
  figInterpreter = defopts(settings, 'Interpreter', 'latex');
  plotDims = defopts(settings, 'PlotDims', dims);
  settings.LegendOption = defopts(settings, 'LegendOption', 'show');
  
  % default function groups
  fcnGroups = {1:5, 6:9, 10:14, 15:19, 20:24, 1:24, ... noiseless
               101:106, 107:121, 122:130, 101:130};   % noisy
  groupNames = {'separable fcns (f1-5)', ... noiseless
                'moderate fcns (f6-9)', ...
                'ill-conditioned fcns (f10-14)', ...
                'multi-modal fcns (f15-19)', ...
                'weakly structured multi-modal fcns (f20-24)', ...
                'all noiseless functions (f1-24)', ...
                'fcns with moderate noise (f101-106)', ... noisy
                'fcns with severe noise (f107-121)', ...
                'highly multi-modal fcns with severe noise (f122-130)', ...
                'all noisy functions (f101-130)'};
  % create groups according to functions to plot
  groupId = cellfun(@(x) any(ismember(x, funs)), fcnGroups);
  fcnGroups = fcnGroups(groupId);
  groupNames = groupNames(groupId);
  
  % settings for all figures
  settings.AggregateFuns = true;
  
  handle = {};
  for d = 1:numel(plotDims)
    % settings for individual dimension
    settings.PlotDims = plotDims(d);
    for fg = 1:numel(fcnGroups)
      % settings for individual function group
      settings.PlotFuns = fcnGroups{fg};
      settings.TitleString = groupNames{fg};
      % remove legendOption if using 'split' option
      removeSplitOption = strcmp(settings.LegendOption, 'split');
      if removeSplitOption
        settings.LegendOption = 'hide';
      end
      % remove legendOption if using 'out' option
      removeOutOption = (strcmp(settings.LegendOption, 'out') && fg < numel(fcnGroups));
      if removeOutOption
        settings.LegendOption = 'hide';
      end
      % plot group
      hdl = relativeFValuesPlot(data, settings);
      % add legend if 'split' option was selected
      if removeSplitOption
        % first plot
        splitEndId = ceil(numOfData/2);
        if isempty(handle)
          % identify marker handles in figure/axes
          legId = find(arrayfun(@(x) ...
            ~strcmp(get(hdl{1}.Children.Children(x), 'Marker'), 'none'), ...
            1:numel(hdl{1}.Children.Children)));
          % first half of names to legend
          legend(hdl{1}.Children.Children(legId(1:splitEndId)), datanames(1:splitEndId), ...
                 'Interpreter', figInterpreter)
        elseif numel(handle) == 1
          % second half of names to legend
          legend(hdl{1}.Children.Children(legId(splitEndId+1:end)), datanames(splitEndId+1:end), ...
                 'Interpreter', figInterpreter)
        end
        settings.LegendOption = 'split';
      end
      % return legendOption if using 'out' option
      if removeOutOption
        settings.LegendOption = 'out';
      end
      % TODO: why not "handle(end+1 : end+numel(hdl)) = hdl;" ???
      if length(hdl) == 2
        handle(end+1:end+2) = hdl;
      else
        handle(end+1) = hdl;
      end
    end
  end
  
end