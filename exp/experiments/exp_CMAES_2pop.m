exp_id = 'exp_CMAES_2pop';
exp_description = 'Pure IPOP-CMA-ES with doubled population size';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 2, 3, 5, 10, 20 }, ...
  'functions',          { 1:24 }, ...           % all functions: num2cell(1:24)
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          { [1:5, 41:50] }, ...   % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }, ...
  'resume',             { true }, ...
};


% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',         { 'none' }, ...         % 'none', 'individual', 'generation', 'restricted'
  'observers',          { [] }, ...             % logging observers
};

modelParams = {};

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/bajeluk/public';
