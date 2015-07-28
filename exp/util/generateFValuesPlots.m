function generateFValuesPlots()
  % function for making graphs showing the dependence of minimal function
  % values on the number of function values
  % Created for GECCO 2015 poster
  
  close all

  % path settings
  exppath = fullfile('exp','experiments');
  gppath = fullfile(exppath,'exp_geneEC_06');
%   rfpath = {[gppath,'_rflite'],[gppath,'_rf5_2'],[gppath,'_rf5_3'],[gppath,'_rf5_4']};
  rfpath = {[gppath,'_rf5_2'],[gppath,'_rf5_3'],[gppath,'_rf5_4']};
  transPath10D = fullfile(exppath,'exp_geneEC_08_10D');
  transPath20D = fullfile(exppath,'exp_geneEC_08_20D');
  cmaespath = fullfile(gppath,'cmaes_results');
  plotResultsFolder = ['/tmp'];
  plotResultsFolder2 = fullfile('..','latex_scmaes','gecco2015poster','images');

  % needed function and dimension settings
  funcSet.BBfunc = [1,2,3,5,6,8,10,11,12,13,14,20,21];
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dims = [2,5,10];
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  
  % loading data
  [trans_evals, trans_settings] = dataready(transPath10D, funcSet, 4);
  [rf_evals, rf_settings] = dataready(rfpath, funcSet, 1, 'rflite');
  [gp_evals, gp_settings] = dataready(gppath, funcSet, 6);
  cmaes_evals = dataready(cmaespath, funcSet, 1, 'cmaes');
  
  % finding data indexes
  set.modelType = 'rf';
  set.evoControlModelGenerations = 1;
  rf1TransId = getSettingsIndex(trans_settings,set);
  rf1Id = getSettingsIndex(rf_settings,set);
  
  set.evoControlModelGenerations = 5;
  rf5TransId = getSettingsIndex(trans_settings,set);
  
  set.modelType = 'gp';
  set.evoControlSampleRange = 1;
  gp5TransId = getSettingsIndex(trans_settings,set);
  
  set.evoControlModelGenerations = 1;
  gp1TransId = getSettingsIndex(trans_settings,set);
  gp1Id = getSettingsIndex(gp_settings,set);
  set.evoControlModelGenerations = 3;
  gp3Id = getSettingsIndex(gp_settings,set);
  
  % color settings
  CMAESCol = [22 22 138];
  GP1TransCol = [255 0 0];
  GP5TransCol = [255 215 0];
  RF1TransCol = [208 32 144];
  RF5TransCol = [0 0 0];
  GP3Col = [100 149 237];
  RF1Col = [116 172 66];
  
%   data = {trans_evals(:,:,rf1TransId),trans_evals(:,:,gp1TransId),trans_evals(:,:,rf5TransId),trans_evals(:,:,gp5TransId),cmaes_evals};
%   datanames = {'RF1T','GP1T','RF5T','GP5T','CMA-ES'};

  data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId),trans_evals(:,:,gp5TransId)};
  datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};
  
%   data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId)};
%   datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};

  
  colors = [CMAESCol; GP3Col; GP5TransCol; RF1Col; RF1TransCol]/255;
  for i = 1:length(funcSet.BBfunc)
    pdfNames{i} = fullfile(plotResultsFolder2,['f',num2str(funcSet.BBfunc(i))]);
  end
  
  han = drawFValuesGraph(data,datanames,funcSet,funcSet.dims,funcSet.BBfunc,colors);
%   print2pdf(han,pdfNames,1)

%   drawGraph(gp_evals,cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [1,2,3,5,6,8];
%   h1 = drawComparison(trans_evals(:,:,1),trans_evals(:,:,2),cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [10,11,12,13,14,20,21];
%   h2 = drawComparison(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);
  
%   pdfname = fullfile(plotResultsFolder,'speedUp');
  print2pdf(han,pdfNames,1)

end

function [data,settings] = dataready(datapath,funcSet,numOfSettings,varargin)
% Prepare data for further processing
% Returns cell array functions x dimensions x settings and appropriate
% settings.
% 
% numOfSettings - number of different settings included in folder
% varargin - directory with multiple s-cmaes settings or pure cmaes (write
% 'cmaes')
  data = cell(length(funcSet.BBfunc),length(funcSet.dims),numOfSettings);
  
  % load and complete results
  if iscell(datapath) % data divided in multiple files
    datalist = {};
    for i = 1:length(datapath)
      list = dir(fullfile(datapath{i},'*.mat'));
      matId = true(1,length(list)); % ids of usable .mat files
      if strfind([list.name],'scmaes_params')
        matId(end) = false;
      end
      if strfind([list.name],'metajob')
        matId(1) = false;
      end    
      datalist(end+1:end+sum(matId)) = {list(matId).name};        % get rid of scmaes_params.mat or metajob.mat
    end
    
    % load data
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
    
    % load settings
    settings = cell(1,numOfSettings);
    for i = 1:numOfSettings
      S = load(datalist{i},'-mat','surrogateParams');
      idx = strfind(datalist{i},'_');
      id = str2double(datalist{i}(1,idx(end)+1:end-4));
      settings{mod(id,numOfSettings)+1} = S.surrogateParams;
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
    
    % load data
    for i = 1:length(datalist)
      S = load(datalist{i},'-mat','y_evals');
%       if ~isempty(varargin) && strcmp(varargin{1},'cmaes')
%         S.y_evals = S.y_evals(1:20); % cmaes ran too many times
%       end
      idx = strfind(datalist{i},'_');
      func = str2double(datalist{i}(1,idx(end-2)+1:idx(end-1)-1)); % function number
      dim = str2double(datalist{i}(1,idx(end-1)+1:idx(end)-2));    % dimension number
      if any(func == funcSet.BBfunc) && any(dim == funcSet.dims)
        id = str2double(datalist{i}(1,idx(end)+1:end-4));            % experiment setting id
        data{funcSet.BBfuncInv(func),funcSet.dimsInv(dim),mod(id,numOfSettings)+1} = S.y_evals;
      end
    end
    
    % load settings
    settings = cell(1,numOfSettings);
    for i = 1:numOfSettings
      S = load(datalist{i},'-mat','surrogateParams');
      idx = strfind(datalist{i},'_');
      id = str2double(datalist{i}(1,idx(end)+1:end-4));
      settings{mod(id,numOfSettings)+1} = S.surrogateParams;
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

function handle = drawFValuesGraph(data,datanames,funcSet,dims,BBfunc,colors)
% Draw graphs of dependence of minimal function values on function 
% evaluations / dimension for individual functions
% data - cell array of data
% datanames - cell array of data names
% funcSet - structure of function and dimension settings
% dims - chosen dimensions
% BBfunc - chosen functions

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
dimIds = funcSet.dimsInv(dims);
funcIds = funcSet.BBfuncInv(BBfunc);

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

function means = gainMeans(data,dimId,funcId,nInstances)
% Returns cell array of means accross chosen dimensions for each function

% cat dimensions if necessary
Dims = length(dimId);
funcs = length(funcId);
means = cell(funcs,1);
if Dims > 1
  for f = 1:funcs
    funcData = [];
    for d = 1:Dims
      actualData = data{funcId(f),dimId(d)};
      useInstances = min([nInstances,size(actualData,2)]);
      funcData = [funcData,actualData(:,1:useInstances)];
    end
    means{f} = mean(funcData,2);
  end
else
  for f = 1:funcs
    actualData = data{funcId(f),dimId};
    useInstances = min([nInstances,size(actualData,2)]);
    means{f} = mean(actualData(:,1:useInstances),2);
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

function index = getSettingsIndex(origSettings, searchSettings)
% returns indeces of searched settings in cell array origSettings
  index = [];
  searchedFields = fieldnames(searchSettings);
  nFields = length(searchedFields);
  searchedValues = cell(1,nFields);
  for j = 1:nFields % find all field values
    searchedValues{j} = getfield(searchSettings,searchedFields{j});
  end
  for i = 1:length(origSettings)
    correctFields = true;
    for j = 1:nFields % compare all needed fields
      correctFields = correctFields && all(getfield(origSettings{i},searchedFields{j}) == searchedValues{j});
    end
    if correctFields
      index(end+1) = i;
    end
  end
end