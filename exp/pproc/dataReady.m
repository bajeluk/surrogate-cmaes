function [data, settings] = dataReady(datapath, funcSet)
% Prepares data for further processing.
% [data, settings] = dataReady(datapath, funcSet) returns cell array 'data'
% of size functions x dimensions x settings and appropriate 'settings'.
%
% Input:
%   datapath - path to data | string
%   funcSet  - structure with fields 'BBfunc' (numbers of BBOB functions) 
%              and 'dims' (numbers of dimensions) | structure
%
% Output:
%   data     - aggregated data of size functions x dimensions x settings 
%              | cell array
%   settings - appropriate settings to 'data', empty if 'datapath' does not
%              contain mat-files with results | structure
%
% See Also:
%   bbobDataReady, catEvalSet

  BBfunc = funcSet.BBfunc;
  dims = funcSet.dims;
  BBfuncInv = inverseIndex(BBfunc);
  dimsInv = inverseIndex(dims);
  nFunc = length(BBfunc);
  nDim  = length(dims);
  
  % load and complete results

  % data are maybe divided between multiple folders
  if (~iscell(datapath))
    datapath = {datapath};
  end
  datalist = {};
  for i = 1:length(datapath)
    actualDataList = gainDataList(datapath{i});
    % find *.mat files according to the IDs (last number before '.mat')
    ids = cellfun(@(x) str2num(x(regexp(x, 'D_\d+\.mat$')+2:end-4)), ...
        actualDataList, 'UniformOutput', false);
    % remove not required files (error, tmp, ...)
    actualDataList = actualDataList(~cellfun(@isempty, ids));
    % sort files according to the IDs
    [~, idsId] = sort(cell2mat(ids));
    actualDataList = actualDataList(idsId);
    datalist(end+1 : end+length(actualDataList)) = actualDataList;
  end
  errPathList = cell2mat(cellfun(@(x) [x, ' '], datapath, 'UniformOutput', false));

  settings = {};
  data = cell(nFunc, nDim);
  if isempty(datalist)
    data = bbobDataReady(datapath, funcSet);
    return
  end
  
  % load data
  for i = 1:length(datalist)
    warning('off', 'MATLAB:load:variableNotFound')
    S = load(datalist{i}, '-mat', 'y_evals', 'surrogateParams', 'cmaesParams');
    warning('on', 'MATLAB:load:variableNotFound')
    if all(isfield(S, {'y_evals', 'surrogateParams', 'cmaesParams'}))
      % unify parameters to one settings structure
      actualSettings = S.surrogateParams;
      fCmaesParams = fields(S.cmaesParams);
      valCmaesParams = struct2cell(S.cmaesParams);
      for f = 1:length(fCmaesParams)
        actualSettings.(fCmaesParams{f}) = valCmaesParams{f};
      end
      % load settings
      settingsId = cellfun(@(x) isequal(actualSettings, x), settings);
      if isempty(settingsId) || ~any(settingsId)
        settings{end+1} = actualSettings;
        settingsId = length(settings);
      else
        settingsId = find(settingsId);
      end

      % load y-evals
      idx = strfind(datalist{i},'_');
      func = str2double(datalist{i}(1,idx(end-2)+1:idx(end-1)-1)); % function number
      dim  = str2double(datalist{i}(1,idx(end-1)+1:idx(end)-2));   % dimension number
      % if we found function and dimension we require, save it to data
      if any(func == BBfunc) && any(dim == dims)
        % if the settings has not yet been loaded for this function and
        % dimension
        if size(data, 3) < settingsId || ...
           isempty(data{BBfuncInv(func), dimsInv(dim), settingsId})
          data{BBfuncInv(func), dimsInv(dim), settingsId} = S.y_evals;
        % else add the result to already existing instances
        else
          data{BBfuncInv(func), dimsInv(dim), settingsId} = ...
            [data{BBfuncInv(func), dimsInv(dim), settingsId}; S.y_evals];
        end
      end
    % necessary variables are missing
    else
      fprintf('Variable ''y_evals'', ''surrogateParams'', or ''cmaesParams'' not found in %s.\n', datalist{i})
    end
  end
  
  % fill remainder of data with empty sets
  if ~isempty(settings) && numel(data) < nFunc*nDim*length(settings)
    data{nFunc, nDim, length(settings)} = [];
  end
    
  data = divSmooth(data, funcSet);
  
end

function datalist = gainDataList(datapath)
% finds list of data in defined path excluding matfiles 'scmaes_params.mat' 
% and 'metajob.mat'

  list = dir(fullfile(datapath, '*.mat'));
  if isempty(list)
    datalist = {};
    return
  end
  % ids of usable .mat files
  matId = true(1, length(list));
  matId = matId & ~cellfun(@(x) strcmp(x, 'scmaes_params.mat'), {list.name});
  matId = matId & cellfun(@(x) isempty(strfind(x, 'metajob.mat')), {list.name});
  % get rid of scmaes_params.mat or metajob.mat
  datalist = cellfun(@(x) fullfile(datapath, x), {list(matId).name}, 'UniformOutput', false);  

end
