function group_results_by_param(exp_id, grouping_param)
% Performs a decomposition of experiments into classes of equality
% w.r.t. all parameters except 'grouping_param'.
%
% Moreover, results of experiments in every class are concatenated
% preserving original ordering.
%
% exp_id           -- unique experiment name
% grouping_param   -- name of parameter over which experiment results will
%                     be concatenated

  if ismember(grouping_param, {'dimensions', 'functions', 'maxfunevals'})
    error('Error: grouping by parameter %s is not supported', grouping_param);
  end

  if ~strcmp(grouping_param, 'instances')
    warning('Note: processing of bbob output implemented only for ''instances'' parameter, continue');
  end

  exppath_short = fileparts(mfilename('fullpath'));

  expFile = fullfile(exppath_short, 'experiments', [exp_id '.m']);
  if ~exist(expFile, 'file')
    error('Error: experiment file "%s" does not exist.', expFile);
  end

  gnuplotScript = 'twoAlgsPlotExtended.gpi';
  gnuplotScript = [exppath_short filesep gnuplotScript];

  % initializes exp_id, bbobParams, surrogateParams, modelParams,
  % cmaesParams among other variables
  run(expFile);

  experimentPath = fullfile(exppath_short, 'experiments', exp_id);

  if ~exist(experimentPath, 'file')
    error('Error: path "%s" not found.', experimentPath);
  end

  % bbob constants
  dataPath = [experimentPath filesep 'bbob_output'];
  indexFileTmpl = 'bbobexp_f%d';
  dataDirTmpl = 'data_f%d';

  % initialize parameter definitions (COPY-PASTED from expInit.m)
  bbParamDef = generateParamDefStruct(bbobParams);
  sgParamDef = generateParamDefStruct(surrogateParams);
  cmParamDef = generateParamDefStruct(cmaesParams);

  % Model settings (is historically part of sgParamDef)
  % -- enumerate all the possible combinations
  sgParamDef(end+1).name = 'modelOpts';
  multiStruct = generateStructOpts(modelParams);
  sgParamDef(end).values = convertMultiStructToCell(multiStruct);

  params = [bbParamDef, sgParamDef, cmParamDef];

  experimentPathGr = fullfile(exppath_short, 'experiments', [exp_id '_grouped_by_' grouping_param]);
  if ~exist(experimentPathGr, 'dir')
    mkdir(experimentPathGr);
  end

  if ~exist([experimentPathGr filesep 'cmaes_results'], 'dir')
    mkdir([experimentPathGr filesep 'cmaes_results']);
  end

  dataPathGr = [experimentPathGr filesep 'bbob_output'];
  if ~exist(dataPathGr, 'dir')
    mkdir(dataPathGr);
  end

  getParamNValues = @(paramDef) cell2mat(structMap(paramDef, @(x) length(x.values)));
  nValues = arrayfun(getParamNValues, [bbParamDef, sgParamDef, cmParamDef]);
  nNonBbobValues = prod(nValues((length(bbParamDef)+1):end));
  nCombinations = prod(nValues);

  % find first field name matching grouping parameter name
  groupingIdx = 1;
  for defIdx = 1:length(params)
    paramDef = params(defIdx);
    if strcmp(paramDef.name, grouping_param)
      break;
    end
    groupingIdx = groupingIdx + 1;
  end

  if groupingIdx > length(nValues)
    error('Error: unknown field name "%s".', grouping_param);
  end

  % number of non-bbob parameters after grouping depends on whether
  % the grouping param itself is member of bbob parameters
  if groupingIdx > length(bbParamDef)
    nNonBbobValuesGr = nNonBbobValues / nValues(groupingIdx);
  else
    nNonBbobValuesGr = nNonBbobValues;
  end

  assert(~mod(nCombinations, nValues(groupingIdx)), 'group_results_by_param(): original number of combinations must be divisible by grouping parameter size');
  nValuesGr = nCombinations / nValues(groupingIdx);

  exp_results = {};
  exp_cmaes_results = {};

  % struct arrays to hold grouped settings & results
  resultsGr = initResults(nValuesGr);
  cmResultsGr = initResults(nValuesGr);
  expSettingsGr = struct('dim', cell(1, nValuesGr), ...
    'bbob_function', cell(1, nValuesGr), ...
    'exp_id', cell(1, nValuesGr), ...
    'instances', cell(1, nValuesGr) ...
  );

  bbParamsGr = initParamsGr(bbParamDef, nValuesGr);
  sgParamsGr = initParamsGr(sgParamDef, nValuesGr);
	cmParamsGr = initParamsGr(cmParamDef, nValuesGr);

  % struct array with index file lines grouped into cell arrays
  bbIndexLinesGr = struct('entry', cell(1, nValuesGr), 'comment', cell(1, nValuesGr), 'data', cell(1, nValuesGr));

  % iterate over original results and group them
  for id = 1:nCombinations
    % ===== OUR RESULTS =====
    % map experiment id to its group
    expIdGr = getGroupExpId(id, nValues, groupingIdx);
    [bbParams, sgParams, cmParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

    % append grouping parameter's value to the group's values
    [bbParams1, sgParams1, cmParams1] = deal(bbParams, sgParams, cmParams);
    if groupingIdx > length(fields(bbParams)) + length(fields(sgParams))
      cmParams1.(grouping_param) = [cmParamsGr(expIdGr).(grouping_param) cmParams1.(grouping_param)];
    elseif groupingIdx > length(fields(bbParams))
      sgParams1.(grouping_param) = [sgParamsGr(expIdGr).(grouping_param) sgParams1.(grouping_param)];
    elseif isfield(bbParamsGr, grouping_param)
      bbParams1.(grouping_param) = [bbParamsGr(expIdGr).(grouping_param) bbParams1.(grouping_param)];
    end

    bbParamsGr(expIdGr) = bbParams1;
    sgParamsGr(expIdGr) = sgParams1;
    cmParamsGr(expIdGr) = cmParams1;

    ifun = bbParams.functions;
    dim = bbParams.dimensions;

    % the first three fields cannot be grouped
    expSettingsGr(expIdGr).dim = dim;
    expSettingsGr(expIdGr).bbob_function = ifun;
    expSettingsGr(expIdGr).maxfunevals = bbParams.maxfunevals;
    expSettingsGr(expIdGr).exp_id = exp_id;
    if strcmp(grouping_param, 'instances')
      expSettingsGr(expIdGr).instances = [expSettingsGr(expIdGr).instances bbParams.instances];
    else
      expSettingsGr(expIdGr).instances = bbParams.instances;
    end

    expFileID = [num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
    resultsFile = [experimentPath filesep exp_id '_results_' expFileID];
    if exist([resultsFile '.mat'], 'file')
      load([resultsFile '.mat'], 'exp_results');
      resultsGr{expIdGr} = resultsGr{expIdGr}.combine(exp_results);
    else
      warning('results file %s.mat does not exist, skipping (exp id %s)', resultsFile, num2str(id));
    end

    cmaesId = floor((id-1) / nNonBbobValues) * nNonBbobValues + 1;
    if id == cmaesId
      cmResultsFile = [experimentPath filesep 'cmaes_results' filesep exp_id '_purecmaes_' expFileID];
      if exist([cmResultsFile '.mat'], 'file')
        load([cmResultsFile '.mat'], 'exp_cmaes_results');
        cmResultsGr{expIdGr} = cmResultsGr{expIdGr}.combine(exp_cmaes_results);
      else
        warning('cmaes results file %s.mat does not exist, skipping (exp id %s)', cmResultsFile, num2str(id));
      end
    end % if

    % ===== BBOB OUTPUT =====
    if ~strcmp(grouping_param, 'instances')
      continue;
    end

    expFileIDGr = [num2str(ifun) '_' num2str(dim) 'D_' num2str(expIdGr)];
    expDataPath = [dataPath filesep expFileID];
    expDataPathGr = [dataPathGr filesep expFileIDGr];

    if ~exist(expDataPathGr, 'dir')
      mkdir(expDataPathGr);
    end

    if ~exist([expDataPathGr filesep sprintf(dataDirTmpl, ifun)], 'dir')
      mkdir([expDataPathGr filesep sprintf(dataDirTmpl, ifun)]);
    end

    indexFile = fullfile(dataPath, expFileID, [sprintf(indexFileTmpl, ifun) '.info']);

    if ~exist(indexFile, 'file')
      warning('index file ''%s'' does not exist, skipping', indexFile);
      continue;
    end

    bbIndexLinesGr = processBbobIndexFile(bbIndexLinesGr, indexFile, expIdGr, expDataPath, expDataPathGr);
   end % for

  % ===== WRITE OUTPUT =====
  for id = 1:numel(resultsGr)
    bbParams = bbParamsGr(id);
    surrogateParams = sgParamsGr(id);
    cmaesParams = cmParamsGr(id);

    exp_settings = expSettingsGr(id);

    ifun = exp_settings.bbob_function;
    dim = exp_settings.dim;
    maxfunevals = exp_settings.maxfunevals;

    results = resultsGr{id};
    if ~results.isEmpty()
      exp_results = results.results;
      y_evals = exp_results.y_evals;
      resultsFile = [experimentPathGr filesep exp_id '_results_' num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
      cmaes_out = [];
      save([resultsFile '.mat'], 'exp_id', 'exp_settings', 'exp_results', 'y_evals', 'surrogateParams', 'cmaesParams', 'bbParams', 'cmaes_out');
    else
      warning('no results for group %d', id);
    end

    cmaesId = floor((id-1) / nNonBbobValuesGr) * nNonBbobValuesGr + 1;
    cmaes_results = cmResultsGr{cmaesId};
    if cmaesId == id
      if ~cmaes_results.isEmpty()
        exp_cmaes_results = cmaes_results.results;
        y_evals = exp_cmaes_results.y_evals;
        cmResultsFile = [experimentPathGr filesep 'cmaes_results' filesep exp_id  '_purecmaes_' num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
        save([cmResultsFile '.mat'], 'exp_id', 'exp_settings', 'exp_cmaes_results', 'y_evals', 'surrogateParams', 'cmaesParams');
      else
        warning('no cmaes results for group %d', id);
      end
    end % if

    % ===== GNUPLOT =====
    if ~results.isEmpty() && ~cmaes_results.isEmpty
      % block COPY-PASTED from bbob_test_01.m
      % Save the data for gnuplot
      gnuplotFile = [experimentPathGr filesep exp_id '_gnuplot_' num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
      generateGnuplotDataExtended([gnuplotFile '.dat'], exp_results, exp_cmaes_results, eval(maxfunevals));

      % save gnuplot script
      if (isfield(surrogateParams, 'modelType')); modelType = surrogateParams.modelType;
      else modelType = ''; end
      if (isfield(surrogateParams, 'evoControl')); evoControl = surrogateParams.evoControl;
      else evoControl = ''; end

      gnuplotScriptCommand = ['sed "s#\<DATAFILE\>#' gnuplotFile '.dat#; s#\<OUTPUTFILE\>#' resultsFile '#; s#\<TITLE\>#f' num2str(ifun) ', ' num2str(dim) 'D#; s#\<DATALINETITLE\>#' modelType ' surrogate, ' evoControl ' EC#; s#\<PARAMS1\>#' sprintfStruct(sgParams, 'escape') '#; s#\<PARAMS2\>#' sprintfStruct(exp_settings, 'escape') '#" ' gnuplotScript ' > ' gnuplotFile '.gpi'];
      disp(gnuplotScriptCommand);
      system(gnuplotScriptCommand);
      % call gnuplot
      system(['gnuplot ' gnuplotFile '.gpi']);

      % print out settings into the text-file
      fid = fopen([resultsFile '.txt'], 'w');
      printSettings(fid, exp_settings, exp_results, surrogateParams, cmaesParams);
      fclose(fid);
    end % if

    % ===== WRITE BBOB INDEX FILES =====
    if ~strcmp(grouping_param, 'instances')
      continue;
    end % fi

    expFileIDGr = [num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
    indexFileGr = fullfile(dataPathGr, expFileIDGr, [sprintf(indexFileTmpl, ifun) '.info']);
    index = bbIndexLinesGr(id);
    algN = min(cellfun(@length, {index.entry, index.comment, index.data}));
    if algN > 0
      indexFileGrId = fopen(indexFileGr, 'w');
      if indexFileGrId < 0
        error('Error: could not open index file ''%s'' for writing', indexFileGr);
      end

      for alg_i = 1:algN
        fprintf(indexFileGrId, '%s\n', index.entry{alg_i});
        fprintf(indexFileGrId, '%s\n', index.comment{alg_i});
        fprintf(indexFileGrId, '%s\n', index.data{alg_i});
      end % for
    end % if
  end % for
end % function


function groupExpId = getGroupExpId(expId, nValues, idx)
% map experiment id to an id after grouping by field 'idx'
% nValues -- a vector with parameter sizes

  % decompose experiment id to parameter indices
  paramIV = getParamIndexVector(expId, nValues);
  % get new bases and indices by leaving out the parameter given by idx
  nValues(idx) = 1;
  paramIV(idx) = 1;
  % re-compose experiment id as a multi-base number with reduced bases and
  % indices
  orders = cumprod(nValues, 'reverse');
  orders = [orders(2:end) 1];
  groupExpId = sum((paramIV - 1) .* orders) + 1;
end % function


function resultsArr = initResults(n)
% initialize an array of results, e.g. one record per task
  resultsArr = cell(1, n);
  for i = 1:numel(resultsArr)
    resultsArr(i) = {ExpResults()};
  end
end


function paramsGr = initParamsGr(paramDef, nValuesGr)
  for i = 1:numel(paramDef)
    paramsGr(nValuesGr).(paramDef(i).name) = [];
  end
end


function c = convertMultiStructToCell(st)
  c = cell(1,length(st));
  for i = 1:length(st)
    c{i} = st(i);
  end
end


function catFiles(fileNames, outName)
% concatenate text files in cell array 'fileIds'; append output to 'outId'
  outId = fopen(outName, 'a');
  if outId < 0
    error('Error: could not open file ''%s'' for append', outName);
  end

  for i = 1:numel(fileNames)
    fileText = fileread(fileNames{i});
    fprintf(outId, '%s', fileText);
  end % for

  fclose(outId);
end % function


function outBbIndexLinesGr = processBbobIndexFile(bbIndexLinesGr, indexFile, expIdGr, expDataPath, expDataPathGr)
% parse index file and save entries for every alg into appropriate group
  indexFileId = fopen(indexFile, 'r');

  alg_i = 0;
  while ~feof(indexFileId)
    alg_i = alg_i + 1;

    % each alg should comprise of three lines
    entryLine = deblank(fgetl(indexFileId));
    assert(~feof(indexFileId));

    if numel(bbIndexLinesGr(expIdGr).entry) < alg_i
      % assume the line ends with algId -- replace exp id in it
      entryLine = regexprep(entryLine, '[0-9]+(_[a-zA-Z-_]*|)''$', sprintf('%d$1''', expIdGr));
      bbIndexLinesGr(expIdGr).entry{alg_i} = entryLine;
    end

    commentLine = deblank(fgetl(indexFileId));
    assert(strfind(commentLine, '%') == 1, 'unexpected input in index file');
    assert(~feof(indexFileId));

    if numel(bbIndexLinesGr(expIdGr).comment) < alg_i
      bbIndexLinesGr(expIdGr).comment{alg_i} = commentLine;
    end

    dataLine = deblank(fgetl(indexFileId));
    data = strsplit(dataLine, ', ');

    if numel(bbIndexLinesGr(expIdGr).data) < alg_i
      % put in the first entry with the filename
      bbIndexLinesGr(expIdGr).data{alg_i} = strjoin(data, ', ');
    else
      % append all entries except the first one which is a filename
      currentEntries = bbIndexLinesGr(expIdGr).data{alg_i};
      newEntries = data(2:end);
      bbIndexLinesGr(expIdGr).data{alg_i} = strjoin([currentEntries, newEntries], ', ');
    end

    % concat data files
    dataFile = fullfile(expDataPath, data{1});
    hdataFile = regexprep(dataFile, '\.dat$', '\.tdat');
    rdataFile = regexprep(dataFile, '\.dat$', '\.rdat');

    if ~exist(dataFile, 'file')
      error('Error: indexed bbob data file %s does not exist', dataFile);
    end

    dataFileGr = fullfile(expDataPathGr, data{1});
    catFiles({dataFile}, dataFileGr);

    if exist(hdataFile, 'file')
      hdataFileGr = regexprep(dataFileGr, '\.dat$', '\.tdat');
      catFiles({hdataFile}, hdataFileGr);
    end

    if exist(rdataFile, 'file')
      rdataFileGr = regexprep(dataFileGr, '\.dat$', '\.rdat');
      catFiles({rdataFile}, rdataFileGr);
    end % if
  end % while

  fclose(indexFileId);

  outBbIndexLinesGr = bbIndexLinesGr;
end % function