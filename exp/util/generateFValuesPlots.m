function generateFValuesPlots()
  % function for making graphs showing the dependence of minimal function
  % values on the number of function values
  % Created for GECCO 2015 poster
  
  close all

  exppath = fullfile('exp','experiments');
  gppath = fullfile(exppath,'exp_geneEC_06');
  rfpath = {[gppath,'_rflite'],[gppath,'_rf5_2'],[gppath,'_rf5_3'],[gppath,'_rf5_4']};
  transPath10D = fullfile(exppath,'exp_geneEC_08_10D');
  transPath20D = fullfile(exppath,'exp_geneEC_08_20D');
  cmaespath = fullfile(gppath,'cmaes_results');
  plotResultsFolder = fullfile('doc','gecco2015paper','images');

  funcSet.BBfunc = [1,2,3,5,6,8,10,11,12,13,14,20,21];
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dims = [2,5,10];
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  
  trans_evals = dataready(transPath10D,funcSet,4);
  rf_evals = dataready(rfpath,funcSet,6,'rflite');
  gp_evals = dataready(gppath,funcSet,6);
  cmaes_evals = dataready(cmaespath,funcSet,1,'cmaes');
  
  hRFall = drawGraph(trans_evals(:,:,4),cmaes_evals,'cmaes',funcSet);
%   pdfname = fullfile('speedUpRF');
%   print2pdf(hRFall,pdfname,1)

%   drawGraph(gp_evals,cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [1,2,3,5,6,8];
%   h1 = drawComparison(trans_evals(:,:,1),trans_evals(:,:,2),cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [10,11,12,13,14,20,21];
%   h2 = drawComparison(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);
  
%   pdfname = fullfile(plotResultsFolder,'speedUp');
%   print2pdf([h1 h2],{[pdfname,'A.pdf'],[pdfname,'B.pdf']},1)

end

function data = dataready(datapath,funcSet,numOfSettings,varargin)
% prepare data for further processing
% numOfSettings - number of different settings included in folder
% varargin - directory with multiple s-cmaes settings or pure cmaes (write
% 'cmaes')
  data = cell(length(funcSet.BBfunc),length(funcSet.dims),numOfSettings);
  
  % load and complete results
  if iscell(datapath) % data divided in multiple files
    datalist = {};
    for i = 1:length(datapath)
      list = dir(fullfile(datapath{i},'*.mat'));
      if i>1
        datalist(end+1:end+length(list)-2) = {list(2:end-1).name};
      else
        datalist(end+1:end+length(list)-1) = {list(1:end-1).name};
      end
    end
    for i = 1:length(datalist)
      idx = strfind(datalist{i},'_');
      id = str2double(datalist{i}(1,idx(end)+1:end-4)); % experiment setting id
      if ~any(strfind(datalist{i},varargin{1})) % if folder does not contain all settings
        id = 1; % TODO: automatically find appropriate setting id
      end
      func = str2double(datalist{i}(1,idx(end-2)+1:idx(end-1)-1)); % function number
      dim = str2double(datalist{i}(1,idx(end-1)+1:idx(end)-2));    % dimension number
      if any(func == funcSet.BBfunc) && any(dim == funcSet.dims)
        S = load(datalist{i},'-mat','y_evals');
        data{funcSet.BBfuncInv(func),funcSet.dimsInv(dim),mod(id,numOfSettings)+1}(end+1:end+length(S.y_evals),1) = S.y_evals;
      end
    end
  else % data in one file
    list = dir(fullfile(datapath,'*.mat'));
    matId = true(1,length(list));
    if strfind([list.name],'scmaes_params')
      matId(end) = false;
    end
    if strfind([list.name],'metajob')
      matId(1) = false;
    end    
    datalist = {list(matId).name};        % get rid of scmaes_params.mat or metajob.mat
    
    for i = 1:length(datalist)
      S = load(datalist{i},'-mat','y_evals');
      if ~isempty(varargin) && strcmp(varargin{1},'cmaes')
        S.y_evals = S.y_evals(1:20); % cmaes ran too many times
      end
      idx = strfind(datalist{i},'_');
      func = str2double(datalist{i}(1,idx(end-2)+1:idx(end-1)-1)); % function number
      dim = str2double(datalist{i}(1,idx(end-1)+1:idx(end)-2));    % dimension number
      if any(func == funcSet.BBfunc) && any(dim == funcSet.dims)
        id = str2double(datalist{i}(1,idx(end)+1:end-4));            % experiment setting id
        data{funcSet.BBfuncInv(func),funcSet.dimsInv(dim),mod(id,numOfSettings)+1} = S.y_evals;
      end
    end
  end
  
  data = divSmooth(data,funcSet);
  
end

function data = divSmooth(data,funcSet)
% divide by dimension and make data smoother
  [func,dims,nSettings] = size(data);
  for s = 1:nSettings
    for d = 1:dims
      for f = 1:func
        fInstant = [];
        for i = 1:length(data{f,d,s})
          data{f,d,s}{i}(:,2) = ceil(data{f,d,s}{i}(:,2)/funcSet.dims(d));
          fInstant(:,end+1) = smoothYEvals(data{f,d,s}{i}(:,1:2),250); % use only first two columns - fitness, evaluations
        end
        data{f,d,s} = fInstant;
      end
    end
  end
end

function handle = drawGraph(data,dataref,refname,funcSet,dims,BBfunc)
% draw graphs of dependence of how better is model than reference data
% according to function evaluations / dimension
% dims - chosen dimensions
% BBfunc - chosen functions

if nargin < 6
  BBfunc = funcSet.BBfunc;
  if nargin < 5
    dims = funcSet.dims;
  end
end

dimIds = funcSet.dimsInv(dims);
funcIds = funcSet.BBfuncInv(BBfunc);

if ~all(dimIds)
  fprintf('Wrong dimesion request!')
end
if ~all(funcIds)
  fprintf('Wrong function request!')
end

% count means
data_means = gainMeans(data,dimIds,funcIds);
dataref_means = gainMeans(dataref,dimIds,funcIds);

% plot results
evaldim = 1:length(data_means{1});
handle = figure();
% add reference line
h(1) = semilogy(evaldim,ones(1,length(evaldim)));
ftitle{1} = refname;
hold on
for f = 1:length(funcIds)
  h(f+1) = semilogy(evaldim,(dataref_means{f}(evaldim))./(data_means{f}(evaldim)));
  ftitle{f+1} = ['f',num2str(BBfunc(f))];
end

ylim(gca,[1e-2 1e2])

legend(h,ftitle,'Location','NorthEastOutside')
xlabel('Number of evaluations / D')
ylabel('\Delta f CMA-ES / \Delta f S-CMA-ES')
hold off

end

function handle = drawComparison(data1,data2,dataref,refname,funcSet,dims,BBfunc)
% draw comparison of two models according to function evaluations / dimension
% returns handle
% dims - chosen dimensions
% BBfunc - chosen functions

if nargin < 7
  BBfunc = funcSet.BBfunc;
  if nargin < 6
    dims = funcSet.dims;
  end
end

dimIds = funcSet.dimsInv(dims);
funcIds = funcSet.BBfuncInv(BBfunc);

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
evaldim = 1:length(dataref_means{1});
scrsz = get(groot,'ScreenSize');
handle = figure('Units','centimeters','Position',[1 scrsz(4)/2 13 7.5]);
subplot(1,2,1);
% add reference line
h(1) = semilogy(evaldim,ones(1,length(evaldim)));
ftitle{1} = refname;
hold on
for f = 1:length(funcIds)
  h(f+1) = semilogy(evaldim,dataref_means{f}(evaldim)./data1_means{f}(evaldim));
  ftitle{f+1} = ['f',num2str(BBfunc(f))];
end
xlabel('Number of evaluations / D')
ylabel('\Deltaf CMA-ES / \Deltaf S-CMA-ES')
legend(h(2:4),ftitle(2:4),'Location','northeast')
title('GP')
ax1 = gca;

subplot(1,2,2);
% add reference line
h(1) = semilogy(evaldim,ones(1,length(evaldim)));
ftitle{1} = refname;
hold on
for f = 1:length(funcIds)
  h(f+1) = semilogy(evaldim,dataref_means{f}(evaldim)./data2_means{f}(evaldim));
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

function means = gainMeans(data,dimId,funcId)
% returns cell array of means accross chosen dimensions

% cat dimensions if necessary
Dims = length(dimId);
funcs = length(funcId);
means = cell(funcs,1);
if Dims > 1
  for f = 1:funcs
    funcData = [];
    for d = 1:Dims
      funcData = [funcData,data{funcId(f),dimId(d)}];
    end
    means{f} = mean(funcData,2);
  end
else
  for f = 1:funcs
    means{f} = mean(data{funcId(f),dimId},2);
  end
end

end

function indInv = inverseIndex(index)
  % returns inversed index vector
  if any(index == 0) % from inversed to normal
    indInv = find(index);
  else % from normal to inversed
    a = zeros(1,max(index));
    a(index) = 1;
    c = cumsum(a);
    a(index) = c(index);
    indInv = a;
  end
end