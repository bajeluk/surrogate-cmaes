function [results, settings, resParams] = metaLearn_loadResults(folders, varargin)
% [results, settings, resParams] = metaLearn_loadResults(folders, varargin) 
% load results of metalearning experiments.
%
% Input:
%   folders  - path to data | string or cell-array of strings
%   varargin - pairs of property (string) and value or struct with 
%              properties as fields:
%     'ExpType'        - type of experiment | {'design', 'cma-run'} 
%     'IgnoreSettings' - list of fields in meta-settings to ignore (e.g. in
%                        case of random initialization) | string or 
%                        cell-array of strings
%     'Instances'      - instances used in the experiment | double
%     'OrigExpName'    - name of the original experiment file | string
%     'SaveResults'    - output file containing loaded results | string
%     'ShowOutput'     - print output of data loading to screen | boolean
%
% Output:
%   results   - aggregated results of size functions x dimensions x 
%               settings | cell-array of structures
%   settings  - appropriate settings to results | cell-array of struct
%   resParams - parameters of loaded results (and experiment settings) |
%               struct
%
% See Also:
%   dataReady, catEvalSet

  if nargout > 0
    results = {};
    settings = {};
  end
  if nargin < 1
    help metaLearn_loadResults
    return
  end
  % cell checking
  if ~iscell(folders)
    folders = {folders};
  end
  function_settings = settings2struct(varargin);
  
  % parse settings
  saveResults = defopts(function_settings, 'SaveResults', '');
  showOutput = defopts(function_settings, 'ShowOutput', false);

  % load data from first folder (temporary)
  [results, settings, resParams] = metaLearn_loadOneFolder(folders{1}, function_settings);
  
  % load data from all folders
%   for s = 1:length(folders)
%     [exp_results{s}, settings{s}] = metaLearn_loadOneFolder(folders{s});
%   end
  
  if ~isempty(saveResults)
    filepath = fileparts(saveResults);
    [~, ~] = mkdir(filepath);
    save(saveResults, 'results', 'settings', 'resParams')
    if showOutput
      fprintf('Results saved to %s\n', saveResults)
    end
  end

end

function [results, settings, res_params] = metaLearn_loadOneFolder(exp_folder, fSettings)
% load results from one folder
  
  % parse settings
  expType = defopts(fSettings, 'ExpType', 'cma-run');

  if strcmp(expType, 'design')
    [results, settings, res_params] = metaLearn_loadDesignedExp(exp_folder, fSettings);
  else
    [results, settings, res_params] = metaLearn_loadRunExp(exp_folder, fSettings);
  end
  
  
end

function [results, settings, res_params] = metaLearn_loadRunExp(exp_folder, fSettings)
% load experiment results from one folder gained through CMA-ES run

  % parse settings
  ignoreSettings = defopts(fSettings, 'IgnoreSettings', {''});
  if ~iscell(ignoreSettings)
    ignoreSettings = {ignoreSettings};
  end
  showOutput = defopts(fSettings, 'ShowOutput', false);
  exp_instances = defopts(fSettings, 'Instances', []);
  origExpLoc = defopts(fSettings, 'OrigExpName', [exp_folder, '.m']);
  
  % gain file list
  if showOutput
    fprintf('Searching for MAT-files in folder %s\n', exp_folder)
  end
  listDir = dir(exp_folder);
  % extract features
  mftId = strcmp({listDir.name}, 'metafeatures');
  mftList = searchFile(fullfile(exp_folder, listDir(mftId).name), '*.mat');
  
  % remove position and dataset folder, checkout feature folder
  remId = strcmp({listDir.name}, '.');
  remId = remId | strcmp({listDir.name}, '..');
  % exclude dataset files (due to large size and therefore long loading
  % time)
  remId = remId | strcmp({listDir.name}, 'dataset');
  remId = remId | mftId;
  listDir(remId) = [];
  folderFiles = cell(numel(listDir), 1);
  % search files in individual folders
  for d = 1:numel(listDir)
    % suppose naming finished by dimension denotation *D.mat to ignore
    % *D_missing_data.mat files
    folderFiles{d} = (searchFile(fullfile(exp_folder, listDir(d).name), '*D.mat'))';
  end
  % cat dataset files
  fileList = [folderFiles{:}]';  
  
  % load original experiment settings if possible
  if isfile(origExpLoc)
    origSettings = getOrigExpSettings(origExpLoc);
    origModelOpts = combineFieldValues(origSettings.modelOptions);
  else
    fprintf('Impossible to load %s. Settings id will be artificial\n', origExpLoc)
  end
  
  % init variables
  settings = {};
%   finished = table();
  finished = [];
  results = table();
  nFiles = numel(fileList);
  % load file sequentially
  for f = 1:nFiles
    if showOutput
      fprintf('Loading %s\n', fileList{f})
    else
      warning('off', 'MATLAB:load:variableNotFound')
    end
    try
      S = load(fileList{f}, '-mat', 'dim', 'fun', 'instances', 'ids', ...
                         'modelOptions', 'modelType', 'stats', 'finished');
    catch err
      warning('%s', getReport(err))
      S = struct();
    end
    warning('on', 'MATLAB:load:variableNotFound')
    if all(isfield(S, {'dim', 'fun', 'instances', 'ids', 'modelOptions', 'stats'}))
      % check problematic loaded variables
      if iscell(S.modelOptions) || ~isfield(S, 'modelType')
        [~, actualFName] = fileparts(fileList{f});
        fname_parts = strsplit(actualFName, '_');
        if numel(fname_parts) ~= 4
          warning('Number of file name parts is not 4 in %s. This may cause some malfunction.', fileList{f})
        end
        % not specified model options will be searched through hash
        % calculation and file name analysis
        if iscell(S.modelOptions)
          loadedHashNumbers = cellfun(@modelHash, S.modelOptions, 'UniformOutput', false);
          mOptsId = strcmp(fname_parts{2}, loadedHashNumbers);
          if sum(mOptsId) == 1
            S.modelOptions = S.modelOptions{mOptsId};
          % catch multiple same hashes
          elseif sum(mOptsId) > 1
            warning('Multiple settings has the same hash number in %s.', fileList{f})
            S.modelOptions = S.modelOptions{find(mOptsId, 1)};
          % no appropriate hash
          else
            warning('No match between hash numbers of model options and file name in %s.', fileList{f})
          end
        end
        % search missing model type in file name
        if ~isfield(S, 'modelType')
          if strcmp(fname_parts{1}(end-4:end), 'model')
            S.modelType = fname_parts{1}(1:end-5);
          else
            S.modelType = fname_parts{1};
          end
        end
      end
      
      % unify parameters to one settings structure
      actualSettings = S.modelOptions;
      actualSettings.modelType = S.modelType;
      actualSettings.modelHash = [S.modelType, '_', modelHash(S.modelOptions)];
      
      % ignore chosen settings
      if ~isempty(ignoreSettings{1})
        actualSettings = safermfield(actualSettings, ignoreSettings);
      end
      
      % get names of statistics
      statNames = fieldnames(S.stats);
      [nInst, ~, nGen] = size(S.stats.(statNames{1}));
      % number of actually calculated ids
      nIds = numel(S.ids);
      % checkout number of instances - if it is not defined (can be found
      % in the experiment), try to estimate it from the rest of data
      % nInst - number of instances in results
      % S.instances - actually calculated results
      if nInst > numel(S.instances)
        if nInst == numel(exp_instances)
          calculatedInst = ismember(exp_instances, S.instances);
          S.instances = exp_instances(calculatedInst);
        % not all instance numbers are known from user -> estimate rows
        elseif isempty(exp_instances) && nInst > numel(exp_instances)
           missResRows = all(isnan(S.stats.(statNames{1})), 2);
           if sum(~missResRows) == nInst
             calculatedInst = ~missResRows;
           else
             warning(['Number of instances in stats is greater than ', ...
                     'number of user defined instances in %s. ', ...
                     'Ignoring extra results.'], fileList{f})
            calculatedInst = [true(numel(exp_instances), 1); false(nInst - numel(exp_instances), 1)];
            S.instances = exp_instances;
           end
        elseif nInst > numel(exp_instances)
          warning(['Number of instances in stats is greater than ', ...
                   'number of user defined instances in %s. ', ...
                   'Ignoring extra results.'], fileList{f})
          calculatedInst = [true(numel(exp_instances), 1); false(nInst - numel(exp_instances), 1)];
          S.instances = exp_instances;
        else
          warning(['Number of instances in stats is lower than ', ...
                   'number of user defined instances in %s. ', ...
                   'Ignoring extra instances.'], fileList{f})
          calculatedInst = true(nInst, 1);
          S.instances = exp_instances(calculatedInst);
        end
      elseif nInst < numel(S.instances)
        warning(['Number of results in stats is lower than ', ...
                 'number of calculated instances in %s. ', ...
                 'Ignoring extra calculated instances.'], fileList{f})
        calculatedInst = true(nInst, 1);
        S.instances = S.instances(1:nInst);
      else
        calculatedInst = true(nInst, 1);
      end
      % recalculate number of instances
      nInst = sum(calculatedInst);
      nActualData = nInst * nIds * nGen;
      
      % find settings identifier
      [settingsId, settings] = propertyId(settings, actualSettings);
      
      % create table of base data identifiers
      actualTable = table(...
                          ones(nActualData, 1) * S.fun, ...
                          ones(nActualData, 1) * S.dim, ...
                          repmat(S.instances(:), nIds*nGen, 1), ...
                          repmat(reshape(repmat(S.ids, nInst, 1), [nInst*nIds, 1]), nGen, 1), ...
                          reshape(repmat(1:nGen, nInst*nIds, 1), [nActualData, 1]), ...
                          ones(nActualData, 1) * settingsId, ...
                          'VariableNames', ...
                          {'fun', 'dim', 'inst', 'id', 'gen', 'model'});
      
      % add statistics to table
      for s = 1:numel(statNames)
        statTable = S.stats.(statNames{s});
        actualTable.(statNames{s}) = reshape(statTable(calculatedInst, 1:nIds, :), ...
                                             [nActualData, 1, 1]);
      end  
      
      % save results to proper position
      results = [results; actualTable];
      
      % notice finished variants
      if isfield(S, 'finished')
        % missing data file
        missingDataFile = [fileList{f}(1:end-4), '_missing_data.mat'];
        if isfile(missingDataFile)
          miss = load(missingDataFile);
          for mi = 1:size(miss.missingData)
            % mark missing data as finished
            S.finished{['id_', num2str(miss.missingData(mi, 4))], ...
                       ['inst_', num2str(miss.missingData(mi, 3))]} = true;
          end
        end
        % only if all settings are finished add them to finished list
        if all(all(table2array(S.finished)))
          % compare to original settings without options added due to
          % result loading
          origSettingsId = propertyId(origModelOpts, ...
            rmfield(actualSettings, {'modelType', 'modelHash'}));
          finished = [finished; S.fun, S.dim, origSettingsId];
        end
      end
      
    end
  end

  % save finished variants
  finishedListFile = fullfile(exp_folder, 'finishedList.txt');
  % open file for writing
  FID = fopen(finishedListFile, 'wt');
  fprintf(FID, '# fun dim model\n');
  % write the matrix
  if FID > 0
    fprintf(FID, ' %d %d %d\n', sortrows(finished)');
    fclose(FID);
  end

  % load metafeatures
  nMFiles = numel(mftList);
  mftTable = table();
  % load metafeature names
  if nMFiles > 0
    % TODO: proper metafeature list (even missing metafeatures considered)
    S = load(mftList{1}, '-mat', 'res');
    % list all metafeature names
    mftsNames = strsplit(printStructure(S.res.ft(1), 'Format', 'field'));
    mftsNames = mftsNames(1:3:end-1);
    mftsNames = cellfun(@(x) strrep(x, '.', '_'), mftsNames', 'UniformOutput', false)';
  end
  % load metafeatures
  for m = 1:nMFiles
    [~, actualFile] = fileparts(mftList{m});
    % gain function, dimension, instance, and id from filename
    fileNumbers = sscanf(actualFile, 'data_f%d_%dD_inst%d_id%d_fts.mat');
    if showOutput
      fprintf('Loading metafeatures in %s\n', mftList{m})
    end
    S = load(mftList{m}, '-mat', 'fun', 'dim', 'inst', 'res');
    % TODO: check matching file name and inner parameters
    nGen = size(S.res.values, 2);
    % prepare metafeature values
    metaVal = mat2cell(S.res.values', nGen, ones(1, numel(mftsNames)));
    % create table of base data identifiers
    actualTable = table(...
                        ones(nGen, 1) * S.fun, ...
                        ones(nGen, 1) * S.dim, ...
                        ones(nGen, 1) * S.inst, ...
                        ones(nGen, 1) * fileNumbers(4), ...
                        (1:nGen)', ...
                        metaVal{:}, ...
                        'VariableNames', ...
                        [{'fun', 'dim', 'inst', 'id', 'gen'}, mftsNames]);
    mftTable = [mftTable; actualTable];
  end
  
  % return loaded functions, dimensions, instances, ids and generations
  res_params.functions = unique(results.fun);
  res_params.dimensions = unique(results.dim);
  res_params.instances = unique(results.inst);
  res_params.ids = unique(results.id);
  res_params.gens = unique(results.gen);
  res_params.mfts = mftTable;
  
end

function [results, settings, res_params] = metaLearn_loadDesignedExp(exp_folder, fSettings)
% load experiment using some sort of design to sample points

  % parse settings
  ignoreSettings = defopts(fSettings, 'IgnoreSettings', {''});
  if ~iscell(ignoreSettings)
    ignoreSettings = {ignoreSettings};
  end
  showOutput = defopts(fSettings, 'ShowOutput', false);

  % load parameter file
  paramFile = fullfile(exp_folder, 'metalearn_params.mat');
  if isfile(paramFile)
    params = load(paramFile);
    if showOutput
      fprintf('Experiment setting file loaded\n')
    end
  else
    warning('Experiment setting file is missing (metalearn_params.mat)')
    params = [];
  end
  % gain file list
  if showOutput
    fprintf('Searching for MAT-files in folder %s\n', exp_folder)
  end
  fileList = searchFile(exp_folder, '*.mat');
  % remove param file
  if ~isempty(params)
    paramsId = strcmp(paramFile, fileList);
    fileList(paramsId) = [];
  end
  
  settings = {};
  results = {};
  instances = {};
  functions = [];
  dimensions = [];
  nFiles = numel(fileList);
  % load file sequentially
  for f = 1:nFiles
    if showOutput
      fprintf('Loading %s\n', fileList{f})
    else
      warning('off', 'MATLAB:load:variableNotFound')
    end
    S = load(fileList{f}, '-mat', 'dim', 'func', 'inst_samp', 'modelOpt', 'modelType', 'results', 'cv');
    warning('on', 'MATLAB:load:variableNotFound')
    if all(isfield(S, {'dim', 'func', 'inst_samp', 'modelOpt', 'modelType', 'results', 'cv'}))
      % unify parameters to one settings structure
      actualSettings = S.modelOpt;
      actualSettings.modelType = S.modelType;
      % if isfield(S, 'inst_samp')
      %   actualSettings.inst_samp = S.inst_samp;
      % end
      
      % ignore chosen settings
      if ~isempty(ignoreSettings{1})
        actualSettings = safermfield(actualSettings, ignoreSettings);
      end
      
      % load dimension, function, and settings
      [dimId, dimensions] = propertyId(dimensions, S.dim);
      [funcId, functions] = propertyId(functions, S.func);
      [settingsId, settings] = propertyId(settings, actualSettings);
      
      % save results to proper position
      results{funcId, dimId, settingsId} = S.results;
      % save instances to proper position
      instances{funcId, dimId, settingsId} = S.inst_samp;

    end
  end
  
  % sort according to function and dimension 
  [functions, funSortId] = sort(functions);
  [dimensions, dimSortId] = sort(dimensions);
  results = results(funSortId, dimSortId, :);
  instances = instances(funSortId, dimSortId, :);
  
  % return loaded functions, dimension and experiment settings
  res_params.functions = functions;
  res_params.dimensions = dimensions;
  res_params.metalearn_params = params;
  res_params.instances = instances;
  
end

function [propId, property] = propertyId(property, value)
% Gain property identifier. If the property value is new, add it to
% property list.
  if iscell(property)
    propId = cellfun(@(x) isequal(value, x), property);
    if isempty(propId) || ~any(propId)
      property{end+1} = value;
      propId = length(property);
    else
      propId = find(propId);
    end
  else
    propId = property == value;
    if isempty(propId) || ~any(propId)
      property(end+1) = value;
      propId = length(property);
    else
      propId = find(propId);
    end
  end
end

function origOpts = getOrigExpSettings(origExpLoc)
% get original experiment settings
  run(origExpLoc)
  clear('origExpLoc')
  variables = who();
  for v = 1:numel(variables)
    origOpts.(variables{v}) = eval(variables{v});
  end
end