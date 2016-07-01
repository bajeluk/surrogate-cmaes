function handle = speedUpPlot(data, dataref, refname, funcSet, dims, BBfunc)
% Plots dependences of how better is model than reference data
% according to function evaluations / dimension.
%
% Input:
%   data    - model data
%   dataref - reference data
%   refname - name of reference data (algorithm)
%   funcSet - structure of function and dimension settings of data
%   dims    - dimensions chosen to plot
%   BBfunc  - functions chosen to plot
%
% Output:
%   handle - handles of resulting figures
%
% See Also:
%   speedUpPlotCompare, fValuesPlot, dataReady

  if nargin < 6
    BBfunc = funcSet.BBfunc;
    if nargin < 5
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
  data_means = gainMeans(data, dimIds, funcIds);
  dataref_means = gainMeans(dataref, dimIds, funcIds);

  % plot results
  evaldim = 1:length(data_means{1});
  handle = figure();
  % add reference line
  h(1) = semilogy(evaldim,ones(1, length(evaldim)));
  ftitle{1} = refname;
  hold on
  for f = 1:length(funcIds)
    h(f+1) = semilogy(evaldim, (dataref_means{f}(evaldim))./(data_means{f}(evaldim)));
    ftitle{f+1} = ['f', num2str(BBfunc(f))];
  end

  ylim(gca, [1e-2, 1e2])

  legend(h, ftitle, 'Location', 'NorthEastOutside')
  xlabel('Number of evaluations / D')
  ylabel('\Delta f CMA-ES / \Delta f S-CMA-ES')
  hold off

end