exp_id = 'exp_GaussianEDA_03';
exp_description = 'Gaussian models EDA algorithm in 27 settings in 2--10D, S-CMA-ES-like experiment';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 2, 3, 5, 10 }, ...
  'functions',          { 1:24 }, ... % all functions: num2cell(1:24)
  'opt_function',       { @opt_GaussianEDA }, ...
  'instances',          { [1:5, 41:50] }, ...   % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }, ...
  'resume',             { false }
};

% Surrogate manager parameters

surrogateParams = { 'xyz', { [] } ...
};

% Model parameters

modelParams = {};

% CEDA parameters

cmaesParams = { ...
  'PopSize',    { '25*dim', '50*dim', '100*dim' }, ...
  'GaussianModelType', { 'GaussianUnivModel', 'GaussianFullModel', 'GaussianNetwork' }, ...
  'SelectionRatio', { 0.1, 0.2, 0.4 };
};

logDir = '/storage/plzen1/home/bajeluk/public';
