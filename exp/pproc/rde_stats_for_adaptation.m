if ~exist('EXPID', 'var')
  EXPID = 'exp_DTSadapt_01'; end

% Histogram bins
bins = 0.025:0.05:1;
% Minimal FE/D to consider
minFEperD = 5;
% Minimal original popsize to consider
minOrigPoints = 1; % 'ceil(0.25 * lambda)';
% Functions
functions = [1:24];
fun_chosen = 1:length(functions);
% Dimensions
dimensions = [2,3,5,10,20];
dim_chosen = [1:5];
% Instances
if ~exist('nInstances', 'var')
  nInstances = 3; end
% Saving pictures
plotPic = false;
savePNG = true;

rnkValid = {};
rnkValidColname = 'rdeValid2';
rnkMeasured = {};
rnkMeasuredColname = 'rdeM1_M2WReplace';
% rnkMeasuredColname = 'rde';
tab = {};
qtab = table([], [], [], [], [], 'VariableNames', {'dim', 'fun', 'Q1', 'Q2', 'Q3'});

% Loading errors from modelStatistics.mat file:
offline = load(['exp/experiments/' EXPID '/modelStatistics.mat']);
rAll = offline.resultTableAll;

for idDim = dim_chosen
  dim = dimensions(idDim);
  ids = ((idDim-1)*(length(functions)*nInstances))+ 1:nInstances:(idDim*(length(functions)*nInstances));
  lambda = 8+floor(6*log(dim));

  for idFun = fun_chosen
    f = functions(idFun);

    %{
    % This is for loading errors from bbob_outpt/*.dat files:
    tab{f} = [];
    instanceToTest = 0;
    while (isempty(tab{f}) && instanceToTest < nInstances)
      tab{f} = readScmaesTxtLog(EXPID, f, dim, ids(idFun) + instanceToTest);
      instanceToTest = instanceToTest + 1;
    end
    if (isempty(tab{f}))
      warning(['Function ' num2str(f) ' in ' num2str(dim) 'D has no data.']);
      rnkValid{f, idDim} = [];
      rnkMeasured{f, idDim} = [];
      continue;
    end

    % take only non-NaN values after minFEperD FE/D
    % and with at least 1/4 of population evaluated
    moreThanQuarter = (tab{f}.orEvalsPop >= myeval(minOrigPoints));
    minEvals = min(minFEperD*dim, tab{f}.totEvals(end-2));
    rnkValid{f, idDim} = tab{f}.rnkValid(tab{f}.totEvals > minEvals & moreThanQuarter & ~isnan(tab{f}.rnkMeasured));
    rnkMeasured{f, idDim} = tab{f}.rnkMeasured(tab{f}.totEvals > minEvals & moreThanQuarter & ~isnan(tab{f}.rnkMeasured));
    %}

    % This is for loading errors from modelStatistics.mat file:
    rnkValid{f, idDim} = rAll{rAll.dim == dim & rAll.fun == f, rnkValidColname};
    rnkValid{f, idDim} = rnkValid{f, idDim}(~isnan(rnkValid{f, idDim}));
    rnkMeasured{f, idDim} = rAll{rAll.dim == dim & rAll.fun == f, rnkMeasuredColname};
    rnkMeasured{f, idDim} = rnkMeasured{f, idDim}(~isnan(rnkMeasured{f, idDim}));
  end

  if (~plotPic)
    continue;
  end

  fig1 = figure();
  fig1.Name = [num2str(dim) 'D'];
  fig1.Position(3) = 600;
  fig1.Position(4) = 1200;

  i = 1;
  for idFun = 1:length(functions)
    f = functions(idFun);
    ax = subplot(length(functions),3,3*(idFun-1)+1);
    % histogram
    hist(rnkMeasured{f, idDim}, bins);
    xlabel('RDE (retrained model)');
    if (length(rnkMeasured{f, idDim}) > 0)
      % normalized histogram values
      h{f} = hist(rnkMeasured{f, idDim}, bins);
      norm_h{f} = h{f}/length(rnkMeasured{f, idDim});
      % quartiles
      hold on;
      quarts = quantile(rnkMeasured{f, idDim}, [0.25, 0.5, 0.75]);
      qtab = [qtab; {dim, f, quarts(1), quarts(2), quarts(3)}];
      for q = quarts
        plot([q q], [max(h{f}), 0], 'r-');
      end
    else
      norm_h{f} = [];
    end
    title(sprintf('Q2=%.2f, Q3=%.2f', quarts(2), quarts(3)));

    ax = subplot(length(functions),3,3*(idFun-1)+2);
    hist(rnkValid{f, idDim}, bins);
    xlabel('RDE (validation set)');
    title(['f' num2str(f)]);

    ax = subplot(length(functions),3,3*idFun);
    scatter(rnkMeasured{f, idDim}, rnkValid{f, idDim});
    if (length(rnkMeasured{f, idDim}) > 0)
      correlation = corr(rnkMeasured{f, idDim}, rnkValid{f, idDim});
    else
      correlation = NaN;
    end
    title(sprintf('corr = %.2f', correlation));
    ax.XLim = [0, 1];
    ax.YLim = [0, 1];
    xlabel('RDE (retrained model)');
    ylabel('RDE (validation set)');
  end

  if (savePNG)
    print([EXPID '_' fig1.Name], '-dpng', '-r80');
  end
end

if (plotPic)
  disp(qtab);
  disp('Average quartiles of RDE from functions 2 and 8 in 2--20D');
  mean(table2array(qtab((qtab.fun == 2 | qtab.fun == 8) & qtab.dim < 40, :)))
  disp('Average quartiles of RDE from function 6 in 2--20D');
  mean(table2array(qtab((qtab.fun == 6) & qtab.dim < 40, :)))
end

%% Calculate ordering of functions based on median validation RDE
tabRnkValidOrdering = table();

for idDim = dim_chosen
  dim = dimensions(idDim);
  medRnkValid = cellfun(@median, rnkValid(:, idDim));
  [~, sortedFunctions] = sort(medRnkValid);
  tmpTable = array2table([sortedFunctions, medRnkValid(sortedFunctions)]);
  tmpTable.Properties.VariableNames = { ['D' num2str(dim)], ['RDEvalid_D' num2str(dim)] };
  % tmpTable = array2table([sortedFunctions, medRnkValid(sortedFunctions)]);
  % tmpTable.Properties.VariableNames = { ['D' num2str(dim) '_rank'], ['D' num2str(dim) '_err'] };  
  tmpTable.Properties.RowNames = arrayfun(@(x) { num2str(x) }, functions);
  tabRnkValidOrdering = [tabRnkValidOrdering, tmpTable];
end

% Calculate median RDE from best functions and worst functions acc. to ValidRDE

nBest = 6;
bestFcn = {@(x) quantile(x, 0.5), @(x) quantile(x, 0.75)};
nWorst = 6;
worstFcn = {@(x) quantile(x, 0.5), @(x) quantile(x, 0.75)};
bestFcns = {};
worstFcns = {};
bestRDEthreshold = NaN(length(bestFcn), length(dim_chosen));
worstRDEthreshold = NaN(length(worstFcn), length(dim_chosen));
bestAggregationFcn = @weightedRankAverage;
worstAggregationFcn = @(x) weightedRankAverage(x, 'reverse');

for idDim = dim_chosen
  dim = dimensions(idDim);
  colTab = ['D' num2str(dim)];
  sortedFcns = tabRnkValidOrdering{:, colTab};
  bestFcns{idDim} = sortedFcns(1:nBest);
  worstFcns{idDim} = sortedFcns((end-nWorst+1):end);
  
  for iBestFcn = 1:length(bestFcn)
    tabColName = ['RDEmeasured_D' num2str(dim) '_' num2str(iBestFcn)];
    rnk2Best{iBestFcn, idDim} = cellfun(bestFcn{iBestFcn}, rnkMeasured(bestFcns{idDim}, idDim));
    warning('off');
    tabRnkValidOrdering{1:nBest, tabColName} = rnk2Best{iBestFcn, idDim};
    warning('on');
    if (any(isnan(rnk2Best{iBestFcn, idDim})))
      warning('Best functions'' rnkMeasured error is NaN for some function');
    end
    % Omit f5 (it has almost no error at all...)
    bestRDEthreshold(iBestFcn, idDim) = bestAggregationFcn(rnk2Best{iBestFcn, idDim}(2:end));
  end

  for iWorstFcn = 1:length(worstFcn)
    tabColName = ['RDEmeasured_D' num2str(dim) '_' num2str(iWorstFcn)];
    rnk2Worst{iWorstFcn, idDim} = cellfun(worstFcn{iWorstFcn}, rnkMeasured(worstFcns{idDim}, idDim));
    warning('off');
    tabRnkValidOrdering{(end-nWorst+1):end, tabColName} = rnk2Worst{iWorstFcn, idDim};
    warning('on');
    if (any(isnan(rnk2Worst{iWorstFcn, idDim})))
      warning('Worst functions'' rnkMeasured error is NaN for some function');
    end
    worstRDEthreshold(iWorstFcn, idDim) = worstAggregationFcn(rnk2Worst{iWorstFcn, idDim});
  end
end
disp([bestRDEthreshold; worstRDEthreshold]);

if (~exist('tabRDEThresholds', 'var')), tabRDEThresholds = table(); end
for iBestFcn = 1:length(bestFcn)
  tabVarName = [EXPID '_best_' num2str(iBestFcn)];
  tabRDEThresholds(:, tabVarName) = num2cell(bestRDEthreshold(iBestFcn,:)');
end
for iWorstFcn = 1:length(worstFcn)
  tabVarName = [EXPID '_worst_' num2str(iWorstFcn)];
  tabRDEThresholds(:, tabVarName) = num2cell(worstRDEthreshold(iWorstFcn,:)');
end

% Let's see, that the RDE quartiles are not possible to express in terms of
% sipmle linear regression of dimension
lmBest = {};
for iBestFcn = 1:length(bestFcn)
  lmBest{iBestFcn} = fitlm(log(dimensions), bestRDEthreshold(iBestFcn, :), 'linear');
end
for iWorstFcn = 1:length(worstFcn)
  lmWorst{iWorstFcn} = fitlm(log(dimensions), worstRDEthreshold(iWorstFcn, :), 'linear');
end

%% take RDE quartiles separately for function best and worst functions
%{
for idDim = dim_chosen
  dim = dimensions(idDim);

  Q_best_nongrupped = qtab(qtab.dim < 40 & ismember(qtab.fun, bestFcns{idDim}),:);
  Q_worst_nongrupped = qtab(qtab.dim < 40 & ismember(qtab.fun, worstFcns{idDim}),:);

  % take average RDE quartiles groupped by dimension
  grouppedbyBest = splitapply(@mean, Q_best_nongrupped{:,3:end}, findgroups(Q_best_nongrupped.dim));
  grouppedbyWorst = splitapply(@mean, Q_worst_nongrupped{:,3:end}, findgroups(Q_worst_nongrupped.dim));
  
  Q_best = Q_best_nongrupped;
  Q_best(Q_best.fun ~= bestFcns{idDim}(1), :) = [];
  Q_best{:, 3:end} = grouppedbyBest;
  Q_wosrt = Q_worst_nongrupped;
  Q_worst(Q_worst.fun ~= bestFcns{idDim}(1), :) = [];
  Q_worst{:, 3:end} = grouppedbyWorst;

  disp('Average quartiles of RDE from functions 2 and 8 for respective dimensions');
  disp(Q_best);
  disp('Average quartiles of RDE from function 6 for respective dimensions');
  disp(Q_worst);

  % Let's see, that the RDE quartiles are not possible to express in terms of
  % sipmle linear regression of dimension
  lmBest = fitlm(Q_best, 'Q2~dim')
  lmWorst = fitlm(Q_worst, 'Q2~dim')
end

%}
