function [data, settings] = dataReady(datapath, funcSet, numOfSettings, varargin)
% Prepares data for further processing.
% Returns cell array 'data' functions x dimensions x settings and 
% appropriate 'settings'.
%
% Input:
%   datapath      - path to data | string
%   funcSet       - structure with fields 'BBfunc' (numbers of BBOB
%                   functions) and 'dims' (numbers of dimensions) 
%                   | structure
%   numOfSettings - number of different settings included in folder 
%                   | integer
%   varargin      - directory with multiple s-cmaes settings or pure cmaes 
%                   (write 'cmaes')
%
% Output:
%   data     - aggregated data of size functions x dimensions x settings 
%              | cell array
%   settings - appropriate settings to 'data' | structure

% TODO: automatically recognize numOfSettings - from scmaes_params.mat

  data = cell(length(funcSet.BBfunc), length(funcSet.dims), numOfSettings);
  
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  
  % load and complete results
  if iscell(datapath) % data divided between multiple folders
    datalist = {};
    for i = 1:length(datapath)
      list = dir(fullfile(datapath{i}, '*.mat'));
      matId = true(1,length(list)); % ids of usable .mat files
      if strfind([list.name], 'scmaes_params')
        matId(end) = false;
      end
      if strfind([list.name], 'metajob')
        matId(1) = false;
      end    
      % get rid of scmaes_params.mat or metajob.mat
      datalist(end+1 : end+sum(matId)) = cellfun(@(x) fullfile(datapath{i}, x), ...
        {list(matId).name}, 'UniformOutput', false); 
    end
    
    % load data
    for i = 1:length(datalist)
      idx = strfind(datalist{i}, '_');
      % experiment setting id
      id = str2double(datalist{i}(1, idx(end)+1:end-4));
      % if folder does not contain all settings
      if ~any(strfind(datalist{i}, varargin{1}))
        % TODO: automatically find appropriate setting id - from scmaes_params.mat
        id = 1;
      end
      func = str2double(datalist{i}(1, idx(end-2)+1 : idx(end-1)-1)); % function number
      dim  = str2double(datalist{i}(1, idx(end-1)+1 : idx(end)-2));   % dimension number
      if any(func == funcSet.BBfunc) && any(dim == funcSet.dims)
        S = load(datalist{i}, '-mat', 'y_evals');
        data{funcSet.BBfuncInv(func), funcSet.dimsInv(dim), mod(id, numOfSettings)+1}(end+1:end+length(S.y_evals), 1) = S.y_evals;
      end
    end
    
    % load settings
    settings = cell(1, numOfSettings);
    for i = 1:numOfSettings
      S = load(datalist{i}, '-mat', 'surrogateParams');
      idx = strfind(datalist{i}, '_');
      id = str2double(datalist{i}(1, idx(end)+1:end-4));
      settings{mod(id, numOfSettings)+1} = S.surrogateParams;
    end
    
  else % data in one folder
    list = dir(fullfile(datapath,'*.mat'));
    matId = true(1,length(list));
    if strfind([list.name],'scmaes_params')
      matId(end) = false;
    end
    if strfind([list.name],'metajob')
      matId(1) = false;
    end    
    % get rid of scmaes_params.mat or metajob.mat
    datalist = cellfun(@(x) fullfile(datapath, x), {list(matId).name}, 'UniformOutput', false);        
    
    % load data
    for i = 1:length(datalist)
      S = load(datalist{i}, '-mat', 'y_evals');
%       if ~isempty(varargin) && strcmp(varargin{1},'cmaes')
%         S.y_evals = S.y_evals(1:20); % cmaes ran too many times
%       end
      idx = strfind(datalist{i},'_');
      func = str2double(datalist{i}(1,idx(end-2)+1:idx(end-1)-1)); % function number
      dim = str2double(datalist{i}(1,idx(end-1)+1:idx(end)-2));    % dimension number
      if any(func == funcSet.BBfunc) && any(dim == funcSet.dims)
        id = str2double(datalist{i}(1,idx(end)+1:end-4));          % experiment setting id
        data{funcSet.BBfuncInv(func),funcSet.dimsInv(dim),mod(id,numOfSettings)+1} = S.y_evals;
      end
    end
    
    % load settings
    settings = cell(1, numOfSettings);
    for i = 1:numOfSettings
      S = load(datalist{i}, '-mat', 'surrogateParams');
      idx = strfind(datalist{i}, '_');
      id = str2double(datalist{i}(1, idx(end)+1:end-4));
      settings{mod(id, numOfSettings) + 1} = S.surrogateParams;
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