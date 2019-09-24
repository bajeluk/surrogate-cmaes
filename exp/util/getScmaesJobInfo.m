function jobInfo = getScmaesJobInfo(exp_id, varargin)
% GETSCMAESJOBINFO - get information about particular jobs from given
% experiment based on scmaes_params.mat file
%
% jobInfo = getScmaesJobInfo(exp_id, jobId) - returns info about given job 
%   number in form of a structure
%
% jobInfo = getScmaesJobInfo(exp_id, 'parameter', value) - returns job ids
% with given value of parameter
%
% jobInfo = ...
% getScmaesJobInfo(exp_id, 'parameter1', value1, 'parameter2', value2, ...) 
%   - returns job ids with given combination of values of appropriate 
%     parameters
%
% Input:
%   exp_id    - experiment identifier (name) | string
%   jobId     - job identifier (1, 2, ..., # setting combinations) |
%               positive scalar integer
%   parameter - settings parameter | string
%   value     - value of appropriate parameter | according to parameter
%
% Output:
%   jobInfo - settings of given job or list of jobIds with given 
%             combination of parameters and values | struct or positive 
%             integer
%
% See Also:
%  getParamIndexVector, getParamsFromIndex

  jobInfo = [];

  if nargin < 1
    help getScmaesJobInfo
    return 
  end
  
  exppath = fullfile('exp', 'experiments');
  paramFile = fullfile(exppath, exp_id, 'scmaes_params.mat');
  assert(isfile(paramFile))

  % load scmaes_params.mat
  params = load(paramFile);
  
  % given job case
  if ~isempty(varargin) && isnumeric(varargin{1})
    jobId = varargin{1};
    [bbParams, sgParams, cmParams] = getParamsFromIndex(jobId, ...
      params.bbParamDef, params.sgParamDef, params.cmParamDef);
    jobInfo = catstruct(bbParams, sgParams, cmParams);
  else
    % TODO
    settings = settings2struct(varargin{:});
    settingFields = fieldnames(settings);
    % get values from firts id
    [bbParams, sgParams, cmParams, ~, nCombs] = getParamsFromIndex(1, ...
      params.bbParamDef, params.sgParamDef, params.cmParamDef);    
    jobStruct = catstruct(bbParams, sgParams, cmParams);
    % add job id if there is match between settings and job info
    if cellfun(@(x) ~any(strcmp(difField(settings, jobStruct), x)), settingFields)
      jobInfo = 1;
    end
    % cycle through all ids and get matches
    for id = 2:nCombs      
      % get values from id
      [bbParams, sgParams, cmParams] = getParamsFromIndex(id, ...
        params.bbParamDef, params.sgParamDef, params.cmParamDef);    
      jobStruct = catstruct(bbParams, sgParams, cmParams);
      % add job id if there is match between settings and job info
      if cellfun(@(x) ~any(strcmp(difField(settings, jobStruct), x)), settingFields)
        jobInfo(end+1) = id;
      end
    end
  end

end
