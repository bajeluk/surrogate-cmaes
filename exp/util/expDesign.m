function setId = expDesign(expName, design, varargin)
% expDesign returns setting ids of model-testing experiment to fulfill
% chosen design.
%
% setId = dOptimalExp(expName, nCombs) returns setting ids for experiment
%   expName for nCombs combinations.
%
% Input:
%   expName - name of the experiment | string
%   design  - chosen experimental design:
%               'dopt' - D-optimal design
%               'full' - full-factorial design
%               'lhs'  - latin-hypercube design
%   settings - pairs of property (string) and value or struct with
%              properties as fields:
%                'Iteration' - number of iterations to achieve 'nCombs' |
%                              positive integer scalar | default: 10^3
%                'nCombs'    - required number of combinations for 'lhs'
%                              design | positive integer scalar
%
% Output:
%   setId - setting ids of experiment (considering only model options) to
%           run to fulfill required design | positive integer vector
%
% See Also:
%   candexch

  if nargout > 0
    setId = [];
  end
  if nargin < 2
    if nargin < 1
      help dOptimalExp
      return
    end
    design = 'lhs';
  end

  % parse input
  settings = settings2struct(varargin);
  expLocation = fullfile('exp', 'experiments', [expName, '.m']);
  assert(isfile(expLocation), ...
    'scmaes:expDesign:wrongFile', ...
    '%s is not a name of an experiment in exp/experiments')
  maxIterations = defoptsi(settings, 'iterations', 10^3);
  assert(isnatural(maxIterations), ...
    'scmaes:expDesign:iterNotNatur', ...
    'Numbar of iterations is not natural.')
  nCombs = defoptsi(settings, 'nCombs', 10);
  assert(isnatural(nCombs), ...
    'scmaes:expDesign:combsNotNatur', ...
    'Number of required combinations is not a natural number.')

  % get model options
  modelOptions = getModelOptions(expName);

  sFields = fieldnames(modelOptions);
  % find modelOptions' cell-array fields
  isCellField = cellfun(@(x) iscell(modelOptions.(x)), sFields);
  % get numbers of settings in cell-array fields
  nCellsField = cellfun(@(x) numel(modelOptions.(x)), sFields(isCellField));
  % get ids of categorical cell-array fields
  isCategoricalField = cellfun(@(x) isCategoricalCell(modelOptions.(x)), ...
                               sFields(isCellField));
  % generate full-factorial design in the reversed order similarly to
  % combineFieldValues
  ffDes = fullfact(nCellsField(end:-1:1));
  % get settings ids
  switch lower(design)
    case {'dopt', 'doptimal', 'd-optimal'}
      % get d-optimal design in the reversed order
      setCombs = dOptimalExp(nCellsField(end:-1:1), isCategoricalField(end:-1:1));
    case {'full', 'fullfactorial', 'full-factorial'}
      setCombs = ffDes;
    case {'lhs', 'lhsdesign'}
      % get lhs design in the reversed order
      setCombs = lhsIterExp(nCombs, nCellsField(end:-1:1), maxIterations);
    otherwise
      error('scmaes:expDesign:wrongdesign', ...
        'There is no design named %s', design)
  end

  % return setId
  setId = find(ismember(ffDes, setCombs, 'rows'));
end

function modelOptions = getModelOptions(expLoc)
  % get model options from experiment definition
  run(expLoc)
  assert(isvarname('modelOptions'), 'scmaes:dOptimalExp:mOptMiss', ...
    'Variable ''modelOptions'' not created during the run of %s', expLoc)
end

function val = isCategoricalCell(cellInput)
% return true when cell is strictly categorical

  % for eval purpose
  dim = 2;

  val = true;
  % try to eval strings to numerical values
  if any(cellfun(@ischar, cellInput))
    for i = 1:numel(cellInput)
      if ischar(cellInput{i})
        % some settings can be evaled later in the code, sometimes dim value
        % is used
        try
          cellInput{i} = eval(cellInput{i});
        catch
          % eval failed => string makes cell-array categorical
          return
        end
      end
    end
  end
  % check if all cells are numeric now
  if all(cellfun(@isnumeric, cellInput))
    val = false;
  end

end

function setCombs = dOptimalExp(nCellsField, isCategoricalField)
% dOptimalExp returns setting combinations of model-testing experiment to
% fulfill D-optimal design

  % settings for cordexch function
  cordexchSettings = {...
    'interaction', ... % d-optimal model type
    'categorical', find(isCategoricalField), ...
    'bounds', arrayfun(@(x) 1:x, nCellsField, 'Uni', false) ...
  };
  % if isempty(nCombs)
  %   setCombs = unique(cordexch(numel(nCellsField), prod(nCellsField), cordexchSettings{:}), 'rows');
  % else
  %   setCombs = cordexch(numel(nCellsField), nCombs, cordexchSettings{:});
  % end
  setCombs = cordexch(numel(nCellsField), prod(nCellsField), cordexchSettings{:});
end

function setCombs = lhsExp(nCombs, nCellsField)
% lhsExp returns setting combinations of model-testing experiment according
% to latin-hypercube sample design
  setCombs = ceil(repmat(nCellsField', nCombs, 1).*lhsdesign(nCombs, numel(nCellsField), 'smooth', 'off'));
end

function setCombs = lhsIterExp(nCombs, nCellsField, maxIterations)
% lhsIterExp returns setting combinations of model-testing experiment
% to latin-hypercube sampling trying to achieve requiered number of
% combinations
  nSetCombs = 0;
  setCombs = [];
  iter = 0;
  fprintf('iteration: combinations\n')
  while nSetCombs < nCombs && iter < maxIterations
    iter = iter + 1;
    newSetCombs = lhsExp(nCombs, nCellsField);
    % number of created combinations
    nNewSetCombs = size(unique(newSetCombs, 'rows'), 1);
    if nNewSetCombs > nSetCombs
      setCombs = newSetCombs;
      nSetCombs = nNewSetCombs;
      fprintf('%9d: %7d\n', iter, nSetCombs);
    end
  end
  % warn if the number of combinations is lower than required
  if nSetCombs < nCombs
    warning('scmaes:expDesign:notAchieveRes', ...
      ['Algorithm did not achieve the required number of %d combinations in %d iterations. ', ...
       'Increase the number of iterations sufficiently to solve this issue.'], ...
      nCombs, maxIterations)
  end
end