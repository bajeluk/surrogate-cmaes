function generateReport(expFolder, varargin)
% generateReport(expFolder, settings) generates report of experiments 
% in expFolder.
%
% Input:
%   expFolder - folder or folders containing experiments (i.e. containing
%               scmaes_params.mat file) | string or cell-array of strings
%   settings - pairs of property (string) and value, or struct with 
%              properties as fields:
%
%     'Description' - description of the report | string
%     'Publish'     - resulting format of published report similar to 
%                     function publish (see help publish) | string
%                   - to disable publishing set option to 'off' (default)
%
% See Also:
%   relativeFValuesPlot, publish

%TODO:
%  - generate report for chosen algorithms
%  - rank table
 
  if nargin < 1
    help generateReport
    return
  end
  
  % parse input
  reportSettings = settings2struct(varargin);
  publishOption = defopts(reportSettings, 'Publish', 'off');
  reportDescription = defopts(reportSettings, 'Description', []);
  if ~iscell(expFolder)
    expFolder = {expFolder};
  end
  
  isFolder = cellfun(@isdir, expFolder);
  assert(any(isFolder), 'generateReport:err:nofolder','No input is a folder')
  % TODO: warning which input folders were not found
  paramFile = cellfun(@(x) fullfile(x, 'scmaes_params.mat'), expFolder(isFolder), 'UniformOutput', false);
  existParFile = cellfun(@(x) logical(exist(x, 'file')), paramFile);
  assert(any(existParFile), 'No input folder contain scmaes_params.mat')
  paramFile = paramFile(existParFile);
  % TODO: warning in which input folders was not found scmaes_params.mat
  
  % initialize key variables
  nParamFiles = length(paramFile);
  settings = cell(nParamFiles, 1);
  expName  = cell(nParamFiles, 1);
  BBfunc   = cell(nParamFiles, 1);
  dims     = cell(nParamFiles, 1);
  % load data
  for s = 1 : nParamFiles
    settings{s} = load(paramFile{s});
    if isfield(settings{s}, 'exp_id')
      expName{s} = settings{s}.exp_id;
    else
      fNameParts = strsplit(paramFile{s}, '/');
      expName{s} = fNameParts{end-1};
    end
    settings{s} = getSettings(settings{s});
    BBfunc{s} = cell2mat(settings{s}.bbParamDef.functions);
    dims{s}   = cell2mat(settings{s}.bbParamDef.dimensions);
  end
  BBfunc = unique([BBfunc{:}]);
  dims = unique([dims{:}]);
  
  % open report file
  ppFolder = cellfun(@(x) fullfile(x(1:end - length([filesep, 'scmaes_params.mat'])), 'pproc'), paramFile, 'UniformOutput', false);
  if ~isdir(ppFolder{1})
    mkdir(ppFolder{1})
  end
  if nParamFiles > 1
    reportName = ['exp_', num2str(nParamFiles), 'report_', num2str(hashGen(expName)), '.m'];
  else
    reportName = [expName{1}, '_report.m'];
  end
  reportFile = fullfile(ppFolder{1}, reportName);
  FID = fopen(reportFile, 'w');
  
  % print report
  
  % introduction
  allExpName = [cellfun(@(x) [x, ', '], expName(1:end-1)', 'UniformOutput', false), expName(end)];
  fprintf(FID, '%%%% %s report\n', [allExpName{:}]);
  if ~isempty(reportDescription)
    fprintf(FID, '%% %s\n', reportDescription);
    fprintf(FID, '%% \n');
  end
  fprintf(FID, '%% Report compares the dependences of minimal function\n');
  fprintf(FID, '%% values on the number of function values of different algorithm settings\n');
  fprintf(FID, '%% tested in experiments %s.\n', [allExpName{:}]);
  fprintf(FID, '%% Moreover, algorithm settings are compared\n');
  fprintf(FID, '%% to important algorithms in continuous black-box optimization field \n');
  fprintf(FID, '%% (CMA-ES, BIPOP-s*ACM-ES, SMAC, S-CMA-ES, and DTS-CMA-ES).\n');
  fprintf(FID, '%% \n');
  for s = 1:nParamFiles
    fprintf(FID, '%% *%s:*\n', expName{s});
    if isfield(settings{s}, 'exp_description')
      fprintf(FID, '%% %s\n', settings{s}.exp_description);
    end
    fprintf(FID, '%% \n');
  end
  fprintf(FID, '%% To gain results publish this script.\n');
  fprintf(FID, '%% \n');
  fprintf(FID, '%% Created on %s', datestr(now));
  if nParamFiles == 1
    fprintf(FID, ' in folder %s', ppFolder{1});
  end
  fprintf(FID, '.\n');
  fprintf(FID, '\n');
  
  % data loading
  fprintf(FID, '%%%% Load data\n');
  fprintf(FID, '\n');
  fprintf(FID, 'expFolder = {};\n');
  for s = 1:nParamFiles
    fprintf(FID, 'expFolder{%d} = ''%s'';\n', s, expFolder{s});
  end
  fprintf(FID, 'resFolder = fullfile(expFolder{1}, ''pproc'');\n');
  fprintf(FID, 'if ~isdir(resFolder)\n');
  fprintf(FID, '  mkdir(resFolder)\n');
  fprintf(FID, 'end\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% loading results\n');
  fprintf(FID, 'funcSet.BBfunc = %s;\n', printStructure(BBfunc, FID, 'Format', 'value'));
  fprintf(FID, 'funcSet.dims = %s;\n', printStructure(dims, FID, 'Format', 'value'));
  fprintf(FID, '[expData, expSettings] = catEvalSet(expFolder, funcSet);\n');
  fprintf(FID, 'nSettings = length(expSettings);\n');
  fprintf(FID, 'expData = arrayfun(@(x) expData(:,:,x), 1:nSettings, ''UniformOutput'', false);\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% create algorithm names\n');
  fprintf(FID, 'expAlgNames = arrayfun(@(x) [''ALG'', num2str(x)], 1:nSettings, ''UniformOutput'', false);\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% color settings\n');
  fprintf(FID, 'expCol = getAlgColors(1:nSettings);\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% load algorithms for comparison\n');
  fprintf(FID, 'algMat = fullfile(''exp'', ''pproc'', ''compAlgMat.mat'');\n');
  fprintf(FID, 'if ~exist(algMat, ''file'')\n');
  fprintf(FID, '  try\n');
  fprintf(FID, '    websave(algMat, ''http://artax.karlin.mff.cuni.cz/~bajel3am/scmaes/compAlgMat.mat'');\n');
  fprintf(FID, '  catch\n');
  fprintf(FID, '  end\n');
  fprintf(FID, 'end\n');
  fprintf(FID, 'try\n');
  fprintf(FID, '  alg = load(algMat);\n');
  fprintf(FID, '  compOn = true;\n');
  fprintf(FID, 'catch\n');
  fprintf(FID, '  compOn = false;\n');
  fprintf(FID, 'end\n');
  fprintf(FID, 'algorithms = alg.algorithm;\n');
  fprintf(FID, '\n');
  
  % experiment settings
  fprintf(FID, '%%%% Experiment settings\n');
  fprintf(FID, '%% \n');
  for s = 1:nParamFiles
    fprintf(FID, '%% * *%s*:\n', expName{s});
    fprintf(FID, '%% \n');
    % keep only parameter fields
    fieldsToRemove = {'exp_id', 'exppath_short', 'exp_description', 'logDir'};
    parSettings = rmfield(settings{s}, fieldsToRemove(isfield(settings{s}, fieldsToRemove)));
    printStructure(parSettings, FID, 'StructName', '%%  ')
    fprintf(FID, '%% \n');
  end
  fprintf(FID, '\n');
  % print algorithms differences
  fprintf(FID, '%% print algorithms differences\n');
  fprintf(FID, 'fprintf(''Algorithm settings differences:\\n\\n'')\n');
  fprintf(FID, '[dFields, dValues] = difField(expSettings);\n');
  fprintf(FID, 'if ~isempty(dFields)\n');
  fprintf(FID, '  for s = 1:nSettings\n');
  fprintf(FID, '    fprintf(''  %%s:\\n'', expAlgNames{s})\n');
  fprintf(FID, '    for f = 1:length(dFields)\n');
  fprintf(FID, '      fprintf(''    %%s = %%s;\\n'', dFields{f}, printStructure(dValues{f, s}, 1, ''Format'', ''value''))\n');
  fprintf(FID, '    end\n');
  fprintf(FID, '    fprintf(''\\n'')\n');
  fprintf(FID, '  end\n');
  fprintf(FID, 'end\n');
  fprintf(FID, '\n');

  
  % tested algorithms comparison
  fprintf(FID, '%%%% Tested algorithms comparison\n');
  fprintf(FID, '\n');
  fprintf(FID, 'for f = funcSet.BBfunc\n');
  fprintf(FID, '  %%%% \n');
  fprintf(FID, '  close all\n');
  fprintf(FID, '  \n');
  fprintf(FID, '  fprintf(''Function %%d\\n'', f)\n');
  fprintf(FID, '  han = relativeFValuesPlot(expData, ...\n');
  fprintf(FID, '                            ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                            ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                            ''PlotFuns'', f, ...\n');
  fprintf(FID, '                            ''AggregateDims'', false, ...\n');
  fprintf(FID, '                            ''AggregateFuns'', false, ...\n');
  fprintf(FID, '                            ''DataNames'', expAlgNames, ...\n');
  fprintf(FID, '                            ''Colors'', expCol, ...\n');
  fprintf(FID, '                            ''LegendOption'', ''out'', ...\n');
  fprintf(FID, '                            ''Statistic'', @median);\n');
  fprintf(FID, 'end\n');
  fprintf(FID, '\n');

  % all algorithms comparison
  fprintf(FID, '%%%% All algorithms comparison\n');
  fprintf(FID, '\n');
  fprintf(FID, 'if compOn\n');
  fprintf(FID, '  data = [expData, {algorithms.data}];\n');
  fprintf(FID, '  datanames = [expAlgNames, {algorithms.name}];\n');
  fprintf(FID, '  colors = [expCol; cell2mat({algorithms.color}'')];\n');
  fprintf(FID, '  \n');
  fprintf(FID, '  for f = funcSet.BBfunc\n');
  fprintf(FID, '    %%%% \n');
  fprintf(FID, '    close all\n');
  fprintf(FID, '    \n');
  fprintf(FID, '    fprintf(''Function %%d\\n'', f)\n');
  fprintf(FID, '    han = relativeFValuesPlot(data, ...\n');
  fprintf(FID, '                              ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                              ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                              ''PlotFuns'', f, ...\n');
  fprintf(FID, '                              ''AggregateDims'', false, ...\n');
  fprintf(FID, '                              ''AggregateFuns'', false, ...\n');
  fprintf(FID, '                              ''DataNames'', datanames, ...\n');
  fprintf(FID, '                              ''Colors'', colors, ...\n');
  fprintf(FID, '                              ''LegendOption'', ''out'', ...\n');
  fprintf(FID, '                              ''Statistic'', @median);\n');
  fprintf(FID, '  end\n');
  fprintf(FID, 'else\n');
  fprintf(FID, '  warning(''Could not load nor download %%s. Omitting comparison of all algorithms.'', algMat);\n');
  fprintf(FID, 'end\n');
  
  % finalize
  fprintf(FID, '\n');
  fprintf(FID, '%%%% Final clearing\n');
  fprintf(FID, 'close all\n');
  fprintf(FID, 'clear s f\n');
  
  % close report file
  fclose(FID);
  
  % copy report file to all pproc folders
  if nParamFiles > 1
    for f = 2 : nParamFiles
      if ~isdir(ppFolder{f})
        mkdir(ppFolder{f})
      end
      copyfile(reportFile, fullfile(ppFolder{f}, reportName));
    end
  end

  % publish report file
  if ~strcmpi(publishOption, 'off')
    fprintf('Publishing %s\nThis may take a few minutes...\n', reportFile)
    addpath(ppFolder{1})
    publishedReport = publish(reportFile, publishOption);
    fprintf('Report published to %s\n', publishedReport)
  end

end

function settings = getSettings(params)
% parses structure fields consisting 'name' and 'value' fields to structure
  paramFields = fieldnames(params);
  for i = 1:length(paramFields)
    if isstruct(params.(paramFields{i})) && all(isfield(params.(paramFields{i}), {'name', 'values'}))
      for f = 1:length(params.(paramFields{i}))
        nnParam = 0;
        if isempty(params.(paramFields{i})(f).name)
          nnParam = nnParam + 1;
          settings.(paramFields{i}).(['NONAMEPARAM', num2str(nnParam)]) = ...
            params.(paramFields{i})(f).values;
        else        
          settings.(paramFields{i}).(params.(paramFields{i})(f).name) = ...
            params.(paramFields{i})(f).values;
        end
      end
    else
      settings.(paramFields{i}) = params.(paramFields{i});
    end
  end
end

function value = getParam(setting, paramName)
% finds and returns appropriate value
  paramID = strcmp({setting.name}, paramName);
  if any(paramID)
    value = setting(paramID).values;
  else
    value = [];
  end
end

function fileNum = hashGen(folders)
% generates hash for result file
 fileNum = num2hex(sum(cellfun(@(x) sum(single(x).*(1:length(x))), folders)));
end