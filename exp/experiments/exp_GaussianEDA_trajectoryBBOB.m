exp_id = 'exp_GaussianEDA_trajectoryBBOB';
exp_description = 'Gaussian EDA algorithm in 18 settings for SAGAS trajectory problem';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 12 }, ...
  'functions',          { 26 }, ...           % all functions: num2cell(1:24)
  'opt_function',       { @opt_GaussianEDA }, ...
  'instances',          { [1:5, 41:50] }, ...   % default is [1:5, 41:50]
  'maxfunevals',        { '30000' }, ...
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
  'PopSize',    { '25*dim', '50*dim', '100*dim' }, ...
  'GaussianModelType', { 'GaussianUnivModel', 'GaussianFullModel', 'GaussianNetwork' }, ...
  'SelectionRatio', { 0.1, 0.2 };
};

logDir = '/storage/plzen1/home/bajeluk/public';
