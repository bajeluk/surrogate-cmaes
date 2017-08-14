exp_id = 'exp_bobyqa_01';
exp_description = 'BOBYQA, S-CMA-ES-like experiment';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 2, 5 }, ...
  'functions',          num2cell(1:24), ...      % all functions: num2cell(1:24)
  'opt_function',       { @opt_bobyqa }, ...
  'instances',          { [1:5, 41:50] }, ...    % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }
};

% Surrogate manager parameters

surrogateParams = { 'xyz', { [] } ...
};

% Model parameters

modelParams = { 'abc', { 0 } ...
};

% BOBYQA parameters

cmaesParams = { ...
  'rho_beg',            { 3.0 }, ...
  'rho_end',            { 1e-7 }, ...
};

logDir = '/storage/plzen1/home/bajeluk/public';
