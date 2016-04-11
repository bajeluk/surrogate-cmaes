function [data, settings] = dataReady(datapath, funcSet)
% Prepares data for further processing.
% Returns cell array 'data' functions x dimensions x settings and 
% appropriate 'settings'.
%
% Input:
%   datapath      - path to data | string
%   funcSet       - structure with fields 'BBfunc' (numbers of BBOB
%                   functions) and 'dims' (numbers of dimensions) 
%                   | structure
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

  % data divided between multiple folders
  if iscell(datapath) 
    datalist = {};
    for i = 1:length(datapath)
      actualDataList = gainDataList(datapath{i});
      datalist(end+1 : end+length(actualDataList)) = actualDataList; 
    end
  % data in one folder
  else 
    datalist = gainDataList(datapath);
  end
  assert(~isempty(datalist), 'Useful data not found')

  settings = {};
  % load data
  for i = 1:length(datalist)
    warning('off', 'MATLAB:load:variableNotFound')
    S = load(datalist{i}, '-mat', 'y_evals', 'surrogateParams');
    warning('on', 'MATLAB:load:variableNotFound')
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
    else
      if ~(isempty(regexp(datalist{i}, '_tmp_\d*.mat', 'once')) || ... % temporary mat-files
          isempty(regexp(datalist{i}, '_\d*_ERROR.mat', 'once')) )    % error files
        fprintf('Variable ''y_evals'' or ''surrogateParams'' not found in %s\n', datalist{i})
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