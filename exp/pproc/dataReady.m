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
      actualDataList = gainDataList(datapath{i});
      datalist(end+1 : end+length(actualDataList)) = actualDataList; 
    end
    assert(~isempty(datalist), 'Useful data not found')
    
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
    settings = loadSettings(datalist, numOfSettings);
    
  else % data in one folder
    datalist = gainDataList(datapath);
    assert(~isempty(datalist), 'Useful data not found')
    
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
    settings = loadSettings(datalist, numOfSettings);
  end
  
  data = divSmooth(data, funcSet);
  
end

function settings = loadSettings(datalist, numOfSettings)
% loads settings from datalist
  settings = cell(1, numOfSettings);
  for i = 1:numOfSettings
    S = load(datalist{i}, '-mat', 'surrogateParams');
    idx = strfind(datalist{i}, '_');
    id = str2double(datalist{i}(1, idx(end)+1:end-4));
    settings{mod(id, numOfSettings) + 1} = S.surrogateParams;
  end
end

function datalist = gainDataList(datapath)
% finds list of data in defined path excluding matfiles 'scmaes_params.mat' 
% and 'metajob.mat'

  list = dir(fullfile(datapath, '*.mat'));
  % ids of usable .mat files
  matId = true(1, length(list));
  if strfind([list.name], 'scmaes_params')
    matId(end) = false;
  end
  if strfind([list.name], 'metajob')
    matId(1) = false;
  end    
  % get rid of scmaes_params.mat or metajob.mat
  datalist = cellfun(@(x) fullfile(datapath, x), {list(matId).name}, 'UniformOutput', false);  

end