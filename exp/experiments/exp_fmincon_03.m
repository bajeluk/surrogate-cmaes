exp_id = 'exp_fmincon_03';
exp_description = 'fmincon Interior-point algorithm with re-starting from optimal point or a uniformly sampled point once in 5 restarts, S-CMA-ES-like experiment';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 2, 3, 5, 10, 20 }, ...
  'functions',          num2cell(1:24), ...           % all functions: num2cell(1:24)
  'opt_function',       { @opt_fmincon }, ...
  'instances',          { [1:5, 41:50] }, ...   % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }, ...
  'resume',             { false }
};

% Surrogate manager parameters

surrogateParams = { 'xyz', { [] } ...
};

% Model parameters

modelParams = {};

% BOBYQA parameters

cmaesParams = { ...
  'Algorithm',          { 'interior-point' }, ...
  ... % 'ObjectiveLimit', ftarget,
  'Display',            { 'off' } ...
};

logDir = '/storage/plzen1/home/bajeluk/public';
