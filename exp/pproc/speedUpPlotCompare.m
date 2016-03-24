function handle = speedUpPlotCompare(data1, data2, dataref, refname, funcSet, dims, BBfunc)
% Draws comparison of two models according to function evaluations / 
% dimension
%
% Input:
%   data1   - data for model 1
%   data2   - data for model 2
%   dataref - reference data
%   refname - name of reference data (algorithm)
%   funcSet - structure of function and dimension settings of data
%   dims    - dimensions chosen to plot
%   BBfunc  - functions chosen to plot
%
% Output:
%   handle - handle of resulting figure
%
% See Also:
%   speedUpPlot, fValuesPlot, dataReady

  if nargin < 7
    BBfunc = funcSet.BBfunc;
    if nargin < 6
      dims = funcSet.dims;
    end
  end

  % get function and dimension IDs
  dimInvIds = inverseIndex(funcSet.dims);
  dimIds = dimInvIds(dims);
  funcInvIds = inverseIndex(funcSet.BBfunc);
  funcIds = funcInvIds(BBfunc);

  if ~all(dimIds)
    fprintf('Wrong dimesion request!')
  end
  if ~all(funcIds)
    fprintf('Wrong function request!')
  end

  % count means
  data1_means = gainMeans(data1,dimIds,funcIds);
  data2_means = gainMeans(data2,dimIds,funcIds);
  dataref_means = gainMeans(dataref,dimIds,funcIds);


  % plot results
  lineWidth = 2;
  evaldim = 1:length(dataref_means{1});
  scrsz = get(groot,'ScreenSize'); % in pixels
  pixPerInch = get(groot, 'ScreenPixelsPerInch');
  inchToCm = 2.54; % 1 inch = 2.54 cm
  scrszCm = scrsz/pixPerInch*inchToCm;
  handle = figure('Units', 'centimeters', 'Position', [1 scrszCm(4)/2 14 8]);
  
  subplot(1,2,1);
  % add reference line
  h(1) = semilogy(evaldim,ones(1,length(evaldim)),'LineWidth',lineWidth);
  ftitle{1} = refname;
  hold on
  for f = 1:length(funcIds)
    h(f+1) = semilogy(evaldim,dataref_means{f}(evaldim)./data1_means{f}(evaldim),'LineWidth',lineWidth);
    ftitle{f+1} = ['f',num2str(BBfunc(f))];
  end
  xlabel('Number of evaluations / D')
  ylabel('\Deltaf CMA-ES / \Deltaf S-CMA-ES')
  legend(h(2:4),ftitle(2:4),'Location','northeast')
  title('GP')
  ax1 = gca;

  subplot(1,2,2);
  % add reference line
  h(1) = semilogy(evaldim,ones(1,length(evaldim)),'LineWidth',lineWidth);
  ftitle{1} = refname;
  hold on
  for f = 1:length(funcIds)
    h(f+1) = semilogy(evaldim,dataref_means{f}(evaldim)./data2_means{f}(evaldim),'LineWidth',lineWidth);
    ftitle{f+1} = ['f',num2str(BBfunc(f))];
  end
  ax2 = gca;
  xlabel('Number of evaluations / D')
  legend(h(5:end),ftitle(5:end),'Location','northeast')
  title('RF')

  % set same axis
  % axYLim = [min([ax1.YLim(1),ax2.YLim(1)]),max([ax1.YLim(2),ax2.YLim(2)])];
  axYLim = [1e-3,max([ax1.YLim(2),ax2.YLim(2)])];
  axXLim = [min(evaldim) max(evaldim)];
  ylim(ax1,axYLim);
  ylim(ax2,axYLim);
  xlim(ax1,axXLim);
  xlim(ax2,axXLim);

  hold off

end