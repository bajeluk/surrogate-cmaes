function reportFile = generateReport(expFolder, varargin)
% reportFile = generateReport(expFolder, settings) generates report of 
% experiments in expFolder.
%
% Input:
%   expFolder - folder or folders containing experiments (i.e. containing
%               scmaes_params.mat file) | string or cell-array of strings
%   settings  - pairs of property (string) and value, or struct with 
%               properties as fields:
%
%     'Description'    - description of the report | string
%     'LegendOption'   - legend option of plots from relativeFValues,
%                        recommended settings for generateReport are:
%                          'out'     - legend is in one separate figure
%                          'manyout' - legend is in multiple separated
%                                      figures
%     'OmitEmptyPlots' - omit plots with no data available | boolean
%     'Publish'        - resulting format of published report similar to 
%                        function publish (see help publish) | string
%                      - to disable publishing set option to 'off' 
%                        (default)
%
% Output:
%   reportFile - name of m-file containing report | string
%
% See Also:
%   relativeFValuesPlot, publish

%TODO:
%  - generate report for chosen algorithms
  
  if nargout > 0
    reportFile = [];
  end
  if nargin < 1
    help generateReport
    return
  end
  
  % parse input
  reportSettings = settings2struct(varargin{:});
  publishOption = defopts(reportSettings, 'Publish', 'off');
  legendOption = defopts(reportSettings, 'LegendOption', 'out');
  omitEmptyPlots = defopts(reportSettings, 'OmitEmptyPlots', false);
  reportDescription = defopts(reportSettings, 'Description', []);
  if ~iscell(expFolder)
    expFolder = {expFolder};
  end
  
  isFolder = cellfun(@isdir, expFolder);
  assert(any(isFolder), 'generateReport:err:nofolder', 'No input is a folder')
  % TODO: warning which input folders were not found
  % remove non-existing folders
  expFolder = expFolder(isFolder);
  paramFile = cellfun(@(x) fullfile(x, 'scmaes_params.mat'), expFolder, 'UniformOutput', false);
  existParFile = cellfun(@(x) logical(exist(x, 'file')), paramFile);
  
  % initialize key variables
  nFolders = sum(isFolder);
  settings = cell(nFolders, 1);
  expName  = cell(nFolders, 1);
  BBfunc   = cell(nFolders, 1);
  dims     = cell(nFolders, 1);
  % evaluations and quantiles to show
  showEval     = [25, 50, 100, 200];
  showQuantile = [0.25, 0.5, 0.75];
  % load data
  for f = 1 : nFolders
    % parametrized experiment
    if existParFile(f)
      settings{f} = load(paramFile{f});
      if isfield(settings{f}, 'exp_id')
        expName{f} = settings{f}.exp_id;
      else
        fNameParts = strsplit(paramFile{f}, filesep);
        expName{f} = fNameParts{end-1};
      end
      settings{f} = getSettings(settings{f});
      BBfunc{f} = cell2mat(settings{f}.bbParamDef.functions);
      dims{f}   = cell2mat(settings{f}.bbParamDef.dimensions);
    % raw data
    else
      folderSplit = strsplit(expFolder{f}, filesep);
      expName{f} = folderSplit{end};
      % extract function and dimension number
      % TODO: speed up this
      tdatFiles = searchFile(expFolder{f}, '*.tdat');
      tdatSplit = strfind(tdatFiles, '_');
      BBfunc{f} = arrayfun(@(x) str2double(tdatFiles{x}(1, tdatSplit{x}(end-1)+2:tdatSplit{x}(end)-1)), ...
        1:length(tdatSplit)); % function numbers
      dims{f} = arrayfun(@(x) str2double(tdatFiles{x}(1, tdatSplit{x}(end)+4:end-5)), ...
        1:length(tdatSplit)); % dimension numbers
    end
  end
  BBfunc = unique([BBfunc{:}]);
  dims = unique([dims{:}]);
  
  % create report name
  if nFolders > 1
    reportName = ['exp_', num2str(nFolders), 'report_', num2str(hashGen(expName)), '.m'];
  else
    reportName = [expName{1}, '_report.m'];
    reportName = repForbiddenChar(reportName, '_');
  end
  % report folder for all generated scripts
  defPpFolder = fullfile('exp', 'pproc', 'generated_scripts');
  if ~isdir(defPpFolder)
    mkdir(defPpFolder)
  end
  % report folder for recent script
  mainPpFolder = fullfile(defPpFolder, reportName(1:end-2));
  if ~isdir(mainPpFolder)
    mkdir(mainPpFolder)
  end
  % open report file
  reportFile = fullfile(mainPpFolder, reportName);
  FID = fopen(reportFile, 'w');
  
  % print report
  
  % introduction
  fprintf(FID, '%%%% %s report\n', strjoin(expName, ', '));
  if ~isempty(reportDescription)
    fprintf(FID, '%% %s\n', reportDescription);
    fprintf(FID, '%% \n');
  end
  fprintf(FID, '%% Report compares the dependences of minimal function\n');
  fprintf(FID, '%% values on the number of function values of different algorithm settings\n');
  fprintf(FID, '%% tested in experiments %s.\n', strjoin(expName, ', '));
  fprintf(FID, '%% Moreover, algorithm settings are compared\n');
  fprintf(FID, '%% to important algorithms in continuous black-box optimization field \n');
  fprintf(FID, '%% (CMA-ES, BIPOP-s*ACM-ES, SMAC, S-CMA-ES, lmm-CMA-ES and DTS-CMA-ES).\n');
  fprintf(FID, '%% \n');
  if any(cellfun(@(x) isfield(x, 'exp_description'), settings))
    for f = 1:nFolders
      fprintf(FID, '%% *%s:*\n', expName{f});
      if isfield(settings{f}, 'exp_description')
        fprintf(FID, '%% %s\n', settings{f}.exp_description);
      end
      fprintf(FID, '%% \n');
    end
  end
  fprintf(FID, '%% Created on %s.\n', datestr(now));
  fprintf(FID, '\n');
  
  % data loading
  fprintf(FID, '%%%%\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% Load data\n');
  fprintf(FID, '\n');
  fprintf(FID, 'expFolder = {};\n');
  for f = 1:nFolders
    fprintf(FID, 'expFolder{%d} = ''%s'';\n', f, expFolder{f});
  end
  fprintf(FID, 'reportLocation = fileparts(which(mfilename));\n');
  fprintf(FID, 'expFolID = strcmp(expFolder, reportLocation);\n');
  fprintf(FID, 'if any(expFolID)\n');
  fprintf(FID, '  resFolder = fullfile(reportLocation, ''pproc'');\n');
  fprintf(FID, 'else\n');
  fprintf(FID, '  resFolder = reportLocation;\n');
  fprintf(FID, 'end\n');
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
  fprintf(FID, '%% create or gain algorithm names\n');
  fprintf(FID, 'expAlgNames = cell(1, nSettings);\n');
  fprintf(FID, 'anonymAlg = ~cellfun(@(x) isfield(x, ''algName''), expSettings);\n');
  fprintf(FID, 'expAlgNames(anonymAlg) = arrayfun(@(x) [''ALG'', num2str(x)], 1:sum(anonymAlg), ''UniformOutput'', false);\n');
  fprintf(FID, 'expAlgNames(~anonymAlg) = cellfun(@(x) x.algName, expSettings(~anonymAlg), ''UniformOutput'', false);\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% color settings\n');
  fprintf(FID, 'expCol = getAlgColors(1:nSettings);\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% evaluation settings\n');
  fprintf(FID, 'showEval = %s;\n', printStructure(showEval, FID, 'Format', 'value'));
  fprintf(FID, '%% quantile settings\n');
  fprintf(FID, 'showQuantile = %s;\n', printStructure(showQuantile, FID, 'Format', 'value'));
  fprintf(FID, '\n');
  fprintf(FID, '%% load algorithms for comparison\n');
  fprintf(FID, '[algData, algNames, algColors] = loadCompAlg(fullfile(''exp'', ''pproc'', ''compAlgMat.mat''), funcSet);\n');
  fprintf(FID, '\n');
  
  % experiment settings
  fprintf(FID, '%%%% Experiment settings\n');
  fprintf(FID, '%% \n');
  for f = 1:nFolders
    fprintf(FID, '%% * *%s*:\n', expName{f});
    fprintf(FID, '%% \n');
    % parametrized algorithm
    if existParFile(f)
      % keep only parameter fields
      fieldsToRemove = {'exp_id', 'exppath_short', 'exp_description', 'logDir'};
      parSettings = rmfield(settings{f}, fieldsToRemove(isfield(settings{f}, fieldsToRemove)));
      printStructure(parSettings, FID, 'StructName', '%%  ')
      fprintf(FID, '%% \n');
    end
  end
  fprintf(FID, '\n');
  % print algorithms differences
  fprintf(FID, '%% print algorithms differences\n');
  fprintf(FID, 'if sum(anonymAlg) > 1\n');
  fprintf(FID, '  fprintf(''Algorithm settings differences:\\n\\n'')\n');
  fprintf(FID, '  [dFields, dValues] = difField(expSettings(anonymAlg));\n');
  fprintf(FID, '  if ~isempty(dFields)\n');
  fprintf(FID, '    for s = 1:nSettings\n');
  fprintf(FID, '      fprintf(''  %%s:\\n'', expAlgNames{s})\n');
  fprintf(FID, '      if anonymAlg(s)\n');
  fprintf(FID, '        for f = 1:length(dFields)\n');
  fprintf(FID, '          fprintf(''    %%s = %%s;\\n'', dFields{f}, printStructure(dValues{f, s}, 1, ''Format'', ''value''))\n');
  fprintf(FID, '        end\n');
  fprintf(FID, '        fprintf(''\\n'')\n');
  fprintf(FID, '      end\n');
  fprintf(FID, '    end\n');
  fprintf(FID, '  end\n');
  fprintf(FID, 'end\n');
  fprintf(FID, '\n');

  
  % tested algorithms comparison
  fprintf(FID, '%%%% Tested algorithms comparison\n');
  fprintf(FID, '\n');
  fprintf(FID, 'close all\n');
  % first rank table
  fprintf(FID, 'rankTable = rankingTable(expData, ''Format'', ''figure'', ...\n');
  fprintf(FID, '                                  ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                                  ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                                  ''DataNames'', expAlgNames, ...\n');
  fprintf(FID, '                                  ''Evaluations'', showEval, ...\n');
  fprintf(FID, '                                  ''Ranking'', ''tolerant'');\n');
  fprintf(FID, '\n');
  fprintf(FID, '%%%%\n');
  fprintf(FID, '%%\n');
  fprintf(FID, '%% *Table 1:* Counts of the 1st ranks of all tested algorithms \n');
  fprintf(FID, '%% from %d benchmark functions \n', length(BBfunc));
  fprintf(FID, '%% according to the lowest achieved ${\\Delta_f}^\\textrm{med}$ for different \n');
  fprintf(FID, '%% FE/D = %s \n', printStructure(showEval, FID, 'Format', 'value'));
  fprintf(FID, '%% and dimensions D = %s.\n', printStructure(dims, FID, 'Format', 'value'));
  fprintf(FID, '%% Ties of the 1st ranks are counted for all respective algorithms. \n');
  fprintf(FID, '%% The ties often occure when $\\Delta f_T = 10^{-8}$ is reached.\n');
  fprintf(FID, '\n');
  % sum of ranks table
  fprintf(FID, 'sumRankTable = rankingTable(expData, ''Format'', ''figure'', ...\n');
  fprintf(FID, '                                     ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                                     ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                                     ''DataNames'', expAlgNames, ...\n');
  fprintf(FID, '                                     ''Evaluations'', showEval, ...\n');
  fprintf(FID, '                                     ''Rank'', ''sum'', ...\n');
  fprintf(FID, '                                     ''Ranking'', ''median'');\n');
  fprintf(FID, '\n');
  fprintf(FID, '%%%%\n');
  fprintf(FID, '%%\n');
  fprintf(FID, '%% *Table 2:* Sums of ranks of all tested algorithms \n');
  fprintf(FID, '%% from %d benchmark functions \n', length(BBfunc));
  fprintf(FID, '%% according to the lowest achieved ${\\Delta_f}^\\textrm{med}$ for different \n');
  fprintf(FID, '%% FE/D = %s \n', printStructure(showEval, FID, 'Format', 'value'));
  fprintf(FID, '%% and dimensions D = %s.\n', printStructure(dims, FID, 'Format', 'value'));
  fprintf(FID, '%% Ties of ranks are replaced by median tied rank for all respective algorithms. \n');
  fprintf(FID, '%% The ties often occure when $\\Delta f_T = 10^{-8}$ is reached.\n');
  fprintf(FID, '%% Missing data ranks are substituted by the average rank (# algorithms + 1)/2.\n');
  fprintf(FID, '%% The sums of ranks are multiplied by factor 10 to show all equalities of ranks.\n');
  fprintf(FID, '\n');
  % EFE table
  fprintf(FID, 'efeTab = efeTable(expData, ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                           ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                           ''DataNames'', expAlgNames, ...\n');
  fprintf(FID, '                           ''Quantiles'', showQuantile);\n');
  fprintf(FID, '\n');
  fprintf(FID, '%%%%\n');
  fprintf(FID, '%%\n');
  fprintf(FID, '%% *Table 3:* Sums of ranks of the expected number of function evaluations ($\\textrm{EFE}$)\n');
  fprintf(FID, '%% of all tested algorithms from %d benchmark functions \n', length(BBfunc));
  fprintf(FID, '%% for different quantiles q = %s \n', printStructure(showQuantile, FID, 'Format', 'value'));
  fprintf(FID, '%% and dimensions D = %s.\n', printStructure(dims, FID, 'Format', 'value'));
  fprintf(FID, '%%\n');
  fprintf(FID, '%% $\\textrm{EFE}$ is calculated as follows:\n');
  fprintf(FID, '%%\n');
  fprintf(FID, '%% $$1)\\ f_{\\textrm{min}} < \\Delta f_T \\qquad \\textrm{EFE} = \\#\\textrm{FE/D}\\,,$$\n');
  fprintf(FID, '%%\n');
  fprintf(FID, '%% $$2)\\ f_{\\textrm{min}} \\geq \\Delta f_T \\qquad \\textrm{EFE} = \n');
  fprintf(FID, '%% \\textrm{FE/D}_\\textrm{max} \\cdot (1 + \\frac{1}{9}\\log \\frac{f_{\\textrm{min}}}{\\Delta f_T})\\,,$$\n');
  fprintf(FID, '%%\n');
  fprintf(FID, '%% where $f_{\\textrm{min}}$ is the best reached function value, $\\Delta f_T\n');
  fprintf(FID, '%% = 10^{-8}$ is target value, $\\textrm{FE/D}$ are function evaluations\n');
  fprintf(FID, '%% divided by dimension, and $\\textrm{FE/D}_\\textrm{max} = 250$.\n');
  fprintf(FID, '%% Missing data ranks are substituted by the average rank (# algorithms + 1)/2.\n');
  fprintf(FID, '%%\n');
  % relativeFValuesPlots
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
  fprintf(FID, '                            ''FunctionNames'', true, ...\n');
  fprintf(FID, '                            ''LegendOption'', ''%s'', ...\n', legendOption);
  fprintf(FID, '                            ''OmitEmpty'', %d, ...\n', omitEmptyPlots);
  fprintf(FID, '                            ''Statistic'', @median);\n');
  fprintf(FID, 'end\n');
  fprintf(FID, '\n');

  % all algorithms comparison
  fprintf(FID, '%%%% All algorithms comparison\n');
  fprintf(FID, '\n');
  fprintf(FID, 'if ~isempty(algData)\n');
  fprintf(FID, '  data = [expData, algData];\n');
  fprintf(FID, '  datanames = [expAlgNames, algNames];\n');
  fprintf(FID, '  colors = [expCol; algColors];\n');
  fprintf(FID, '  \n');
  fprintf(FID, '  close all\n');
  % first rank table
  fprintf(FID, '  rankTable = rankingTable(data, ''Format'', ''figure'', ...\n');
  fprintf(FID, '                                 ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                                 ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                                 ''DataNames'', datanames, ...\n');
  fprintf(FID, '                                 ''Evaluations'', showEval, ...\n');
  fprintf(FID, '                                 ''Ranking'', ''tolerant'');\n');
  fprintf(FID, '  \n');
  fprintf(FID, '  %%%%\n');
  fprintf(FID, '  %%\n');
  fprintf(FID, '  %% *Table 4:* Counts of the 1st ranks of all compared algorithms \n');
  fprintf(FID, '  %% from %d benchmark functions \n', length(BBfunc));
  fprintf(FID, '  %% according to the lowest achieved ${\\Delta_f}^\\textrm{med}$ for different \n');
  fprintf(FID, '  %% FE/D = %s \n', printStructure(showEval, FID, 'Format', 'value'));
  fprintf(FID, '  %% and dimensions D = %s.\n', printStructure(dims, FID, 'Format', 'value'));
  fprintf(FID, '  %% Ties of the 1st ranks are counted for all respective algorithms. \n');
  fprintf(FID, '  %% The ties often occure when $\\Delta f_T = 10^{-8}$ is reached.\n');
  fprintf(FID, '  \n');
  % sum of ranks table
  fprintf(FID, '  sumRankTable = rankingTable(data, ''Format'', ''figure'', ...\n');
  fprintf(FID, '                                    ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                                    ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                                    ''DataNames'', datanames, ...\n');
  fprintf(FID, '                                    ''Evaluations'', showEval,...\n');
  fprintf(FID, '                                    ''Rank'', ''sum'', ...\n');
  fprintf(FID, '                                    ''Ranking'', ''median'');\n');
  fprintf(FID, '  \n');
  fprintf(FID, '  %%%%\n');
  fprintf(FID, '  %%\n');
  fprintf(FID, '  %% *Table 5:* Sums of ranks of all compared algorithms \n');
  fprintf(FID, '  %% from %d benchmark functions \n', length(BBfunc));
  fprintf(FID, '  %% according to the lowest achieved ${\\Delta_f}^\\textrm{med}$ for different \n');
  fprintf(FID, '  %% FE/D = %s \n', printStructure(showEval, FID, 'Format', 'value'));
  fprintf(FID, '  %% and dimensions D = %s.\n', printStructure(dims, FID, 'Format', 'value'));
  fprintf(FID, '  %% Ties of ranks are replaced by median tied rank for all respective algorithms. \n');
  fprintf(FID, '  %% The ties often occure when $\\Delta f_T = 10^{-8}$ is reached.\n');
  fprintf(FID, '  %% Missing data ranks are substituted by the average rank (# algorithms + 1)/2.\n');
  fprintf(FID, '  %% The sums of ranks are multiplied by factor 10 to show all equalities of ranks.\n');
  fprintf(FID, '  \n');
  % EFE table
  fprintf(FID, '  efeTab = efeTable(data, ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                          ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                          ''DataNames'', datanames, ...\n');
  fprintf(FID, '                          ''Quantiles'', showQuantile);\n');
  fprintf(FID, '  \n');
  fprintf(FID, '  %%%%\n');
  fprintf(FID, '  %%\n');
  fprintf(FID, '  %% *Table 6:* Sums of ranks of the expected number of function evaluations ($\\textrm{EFE}$)\n');
  fprintf(FID, '  %% of all compared algorithms from %d benchmark functions \n', length(BBfunc));
  fprintf(FID, '  %% for different quantiles q = %s \n', printStructure(showQuantile, FID, 'Format', 'value'));
  fprintf(FID, '  %% and dimensions D = %s.\n', printStructure(dims, FID, 'Format', 'value'));
  fprintf(FID, '  %%\n');
  fprintf(FID, '  %% $\\textrm{EFE}$ is calculated as follows:\n');
  fprintf(FID, '  %%\n');
  fprintf(FID, '  %% $$1)\\ f_{\\textrm{min}} < \\Delta f_T \\qquad \\textrm{EFE} = \\#\\textrm{FE/D}\\,,$$\n');
  fprintf(FID, '  %%\n');
  fprintf(FID, '  %% $$2)\\ f_{\\textrm{min}} \\geq \\Delta f_T \\qquad \\textrm{EFE} = \n');
  fprintf(FID, '  %% \\textrm{FE/D}_\\textrm{max} \\cdot (1 + \\frac{1}{9}\\log \\frac{f_{\\textrm{min}}}{\\Delta f_T})\\,,$$\n');
  fprintf(FID, '  %%\n');
  fprintf(FID, '  %% where $f_{\\textrm{min}}$ is the best reached function value, $\\Delta f_T\n');
  fprintf(FID, '  %% = 10^{-8}$ is target value, $\\textrm{FE/D}$ are function evaluations\n');
  fprintf(FID, '  %% divided by dimension, and $\\textrm{FE/D}_\\textrm{max} = 250$.\n');
  fprintf(FID, '  %% Missing data ranks are substituted by the average rank (# algorithms + 1)/2.\n');
  fprintf(FID, '  %%\n');
  % relativeFValuesPlots
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
  fprintf(FID, '                              ''FunctionNames'', true, ...\n');
  fprintf(FID, '                              ''LegendOption'', ''%s'', ...\n', legendOption);
  fprintf(FID, '                              ''OmitEmpty'', %d, ...\n', omitEmptyPlots);
  fprintf(FID, '                              ''Statistic'', @median);\n');
  fprintf(FID, '  end\n');
  fprintf(FID, '  \n');
  % summary graphs
  fprintf(FID, '  %%%% Summary graphs\n');
  fprintf(FID, '  %% Summary graphs are averaged through all functions for individual algorithms in separate dimensionalities.\n');
  fprintf(FID, '  close all\n');
  fprintf(FID, '  han = relativeFValuesPlot(data, ...\n');
  fprintf(FID, '                          ''DataDims'', funcSet.dims, ...\n');
  fprintf(FID, '                          ''DataFuns'', funcSet.BBfunc, ...\n');
  fprintf(FID, '                          ''AggregateDims'', false, ...\n');
  fprintf(FID, '                          ''AggregateFuns'', true, ...\n');
  fprintf(FID, '                          ''DataNames'', datanames, ...\n');
  fprintf(FID, '                          ''Colors'', colors, ...\n');
  fprintf(FID, '                          ''FunctionNames'', true, ...\n');
  fprintf(FID, '                          ''LegendOption'', ''%s'', ...\n', legendOption);
  fprintf(FID, '                          ''OmitEmpty'', %d, ...\n', omitEmptyPlots);
  fprintf(FID, '                          ''Statistic'', @median);\n');
  fprintf(FID, 'else\n');
  fprintf(FID, '  warning(''Could not load %%s.\\n');
  fprintf(FID, 'For the latest version download http://artax.karlin.mff.cuni.cz/~bajel3am/scmaes/compAlgMat.mat\\n');
  fprintf(FID, 'Omitting comparison of all algorithms.'', algMat);\n');
  fprintf(FID, 'end\n');
  
  % finalize
  fprintf(FID, '\n');
  fprintf(FID, '%%%%\n');
  fprintf(FID, '\n');
  fprintf(FID, '%% Final clearing\n');
  fprintf(FID, 'close all\n');
  fprintf(FID, 'clear s f\n');
  
  % close report file
  fclose(FID);
  fprintf('Report generated to %s\n', reportFile)
  
  % copy report file to all folders of parametrized algorithms
  parFolders = expFolder(existParFile);
  for f = 1 : sum(existParFile)
    ppFolder = fullfile(parFolders{f}, 'pproc');
    if ~isdir(ppFolder)
      mkdir(ppFolder)
    end
    copyfile(reportFile, fullfile(ppFolder, reportName));
  end

  % publish report file
  if ~strcmpi(publishOption, 'off')
    % warn user about clearing variables
    fprintf('Recommendation: Use ')
    fprintf(2, '''clear all'' ');
    fprintf('command before running generateReport to ensure that no static variables from the previous run will be used.\n')
    fprintf('Publishing...\nThis may take a few minutes...\n')
    addpath(mainPpFolder)
    publishedReport = publish(reportFile, 'format', publishOption, ...
                                          'showCode', false);
    % replace importhtml tags with generated html code
    if strcmp(publishOption, 'html')
      htmlPostProc(publishedReport)
    end
    fprintf('Report published to %s\n', publishedReport)
  end
  
  % return report file only if needed
  if nargout < 1
    clear reportFile
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

function str = repForbiddenChar(str, newChar)
% replace forbidden characters by new character
  forbidden = {'-', ' '};
  for f = 1:length(forbidden)
    str = strrep(str, forbidden{f}, newChar);
  end
end
