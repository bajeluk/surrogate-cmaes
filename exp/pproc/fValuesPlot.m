function handle = fValuesPlot(data, datanames, funcSet, dims, BBfunc, colors)
% Plots dependences of minimal function values on function 
% evaluations / dimension for individual functions.
%
% Input:
%   data      - cell array of data
%   datanames - cell array of data names (e.g. names of algorithms)
%   funcSet   - structure of function and dimension settings of data
%   dims      - chosen dimensions
%   BBfunc    - chosen BBOB functions
%   colors    - colors of individual functions
%
% Output:
%   handle - handles of resulting figures
%
% See Also:
%   speedUpPlot, speedUpPlotCompare, dataReady

  numOfData = length(datanames);

  if nargin < 6
    colors = rand(numOfData,3);
    if nargin < 5
      BBfunc = funcSet.BBfunc;
      if nargin < 4
        dims = funcSet.dims;
      end
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
  data_means = cellfun(@(D) gainMeans(D,dimIds,funcIds,useMaxInstances),data,'UniformOutput',false);

  % plot results
  evaldim = 1:length(data_means{1}{1});
  linewidth = 2;

  for f = 1:length(funcIds)
    handle(f) = figure('Units','centimeters','Position',[1 1 12.5 6]);
    h(1) = semilogy(evaldim,data_means{1}{f}(evaldim),'LineWidth',linewidth,'Color',colors(1,:));
    ftitle{1} = datanames{1};
    hold on
    grid on
    for D = 2:length(datanames)
      h(D) = semilogy(evaldim,data_means{D}{f}(evaldim),'LineWidth',linewidth,'Color',colors(D,:));
      ftitle{D} = datanames{D};
    end

    % additional plot settings
  %   ylim(gca,[1e-8 1e5])

    legend(h,ftitle,'Location','NorthEast')
    title(['f',num2str(funcSet.BBfunc(f))])
    xlabel('Number of evaluations / D')
    ylabel('Minimum function values')
    hold off
  end

end