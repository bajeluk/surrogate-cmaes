function nlp = infPriorOnly(inf, prior, hyp, varargin)
% repjak, 2018: Just evaluates the prior.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Hannes Nickisch and Roman Garnett, 2016-10-26.
%
% See also INFMETHODS.M, USAGEPRIOR.M, PRIORDISTRIBUTIONS.M.

nlp = 0;
if ~isempty(prior)                % evaluate hyperprior
  nam = fieldnames(prior);
  num = zeros(numel(any2vec(hyp)),1);      % number of priors per hyperparameter
  for i=1:numel(nam)     % iterate over kinds of hyperparameters cov/lik/mean/xu
    ni = nam{i};                                         % name of the ith field
    if strcmp(ni,'multi')                      % first catch multivariate priors
      p = prior.(ni);
      for j=1:numel(p)  % iterate over individual muti-hyperparameter components
        pj = p{j};
        if ~isempty(pj)          % only proceed if a nonempty prior is specified
          idx = pj{end}; pj = pj(1:end-1);            % grab index from cell end
          pjstr = pj{1}; if ~ischar(pjstr), pjstr = func2str(pjstr); end
          if numel(pjstr)<5 || ~strcmp(pjstr(end-4:end),'Multi')
            error('multivariate hyperpriors are called <Name>Multi')
          end
          if isstruct(idx)                           % massage idx into a vector
            idxs = vec2any(hyp,zeros(size(num)));             % structured index
            for nj = fieldnames(idx), idxs.(nj{1})( idx.(nj{1}) ) = 1; end
            idx = any2vec(idxs)>0;                  % linearise structured index
          else
            idxz = zeros(size(num)); idxz(idx) = 1; idx = idxz>0; % binary index
          end
          if sum(idx)<=1, error('multivariate priors need >1 hyperparams'), end
          num(idx) = num(idx)+1;                             % inc prior counter
          hypu = any2vec(hyp);
          lp = feval(pj{:}, hypu(idx));    % evaluate prior distribution
          nlp = nlp-lp;
        else
          error('multivariate priors should be non empty')
        end
      end
      continue                                       % jump to univariate priors
    end
    if ~isfield(hyp,ni), error(['unknown prior field ',ni,' in hyp']), end
    p = prior.(ni);
    if numel(p)~=numel(hyp.(ni)), error(['bad hyp/prior field length ',ni]), end
    for j=1:numel(p)         % iterate over individual hyperparameter components
      pj = p{j}; if ~iscell(pj) && ~isempty(pj), pj = {pj}; end   % enforce cell
      if ~isempty(pj)            % only proceed if a nonempty prior is specified
        num = vec2any(hyp,num); num.(ni)(j) = num.(ni)(j)+1; % inc prior counter
        num = any2vec(num);
        pj1str = pj{1}; if ~ischar(pj1str), pj1str = func2str(pj1str); end
        lp = feval(pj{:}, hyp.(ni)(j));    % evaluate prior distribution
        nlp = nlp - lp;
      end
    end
  end
  if any(num>1)                 % check for hypers with more than a single prior
    num = vec2any(hyp,num); nam = fieldnames(num);
    s = '';
    for i=1:numel(nam)
      idx = find(num.(nam{i})>1);
      for j=1:numel(idx)
        s = [s,sprintf('hyp.%s(%d) ',nam{i},idx(j))];
      end
    end
    error(['More than 1 prior specified for ',s(1:end-1),'.'])  
  end
end