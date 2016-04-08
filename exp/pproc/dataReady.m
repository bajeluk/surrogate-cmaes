function [data, settings] = dataReady(datapath, funcSet,  varargin)
% Prepares data for further processing.
% Returns cell array 'data' functions x dimensions x settings and 
% appropriate 'settings'.
%
% Input:
%   datapath      - path to data | string
%   funcSet       - structure with fields 'BBfunc' (numbers of BBOB
%                   functions) and 'dims' (numbers of dimensions) 
%                   | structure
%   varargin      - directory with multiple s-cmaes settings or pure cmaes 
%                   (write 'cmaes')
%
% Output:
%   data     - aggregated data of size functions x dimensions x settings 
%              | cell array
%   settings - appropriate settings to 'data' | structure

  BBfunc = funcSet.BBfunc;
  dims = funcSet.dims;
  BBfuncInv = inverseIndex(BBfunc);
  dimsInv = inverseIndex(dims);
  
  data = cell(length(BBfunc), length(dims));
  
  % load and complete results
  if iscell(datapath) % data divided between multiple folders
    datalist = {};
    for i = 1:length(datapath)
      actualDataList = gainDataList(datapath{i});
      datalist(end+1 : end+length(actualDataList)) = actualDataList; 
    end
  else % data in one folder
    datalist = gainDataList(datapath);
  end
  assert(~isempty(datalist), 'Useful data not found')

  settings = {};
  % load data
  for i = 1:length(datalist)
    S = load(datalist{i}, '-mat', 'y_evals', 'surrogateParams');
    if all(isfield(S, {'y_evals', 'surrogateParams'}))
      % load settings
      settingsId = cellfun(@(x) isequal(S.surrogateParams, x), settings);
      if isempty(settingsId) || ~any(settingsId)
        settings{end+1} = S.surrogateParams;
        settingsId = length(settings);
      else
        settingsId = find(settingsId);
      end

      % load y-evals
      idx = strfind(datalist{i},'_');
      func = str2double(datalist{i}(1,idx(end-2)+1:idx(end-1)-1)); % function number
      dim  = str2double(datalist{i}(1,idx(end-1)+1:idx(end)-2));   % dimension number
      if any(func == BBfunc) && any(dim == dims)
        data{BBfuncInv(func), dimsInv(dim), settingsId} = S.y_evals;
      end
    end
  end
    
  data = divSmooth(data, funcSet);
  
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