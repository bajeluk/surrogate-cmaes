% VAL = BENCHMARKSEXTRA(X, FUNCID)
% VAL = BENCHMARKSEXTRA(X, STRFUNC)
%    Input: 
%       X -- solution column vector or matrix of column vectors
%       FUNCID -- number of function to be executed with X as input
%       STRFUNC -- function as string to be executed with X as input
%    Output: function value(s) of solution(s)
%    Examples: 
%      F = BENCHMARKSEXTRA([1 2 3]', 117); 
%      F = BENCHMARKSEXTRA([1 2 3]', 'f101'); 
% 
% NBS = BENCHMARKSEXTRA() 
% NBS = BENCHMARKSEXTRA('FunctionIndices') 
%    Output: 
%      NBS -- array of valid benchmark function numbers, 
%             presumably 101:130
%
% FHS = BENCHMARKSEXTRA('handles')
%    Output: 
%      FHS -- cell array of function handles
%
% STR = BENCHMARKSEXTRA('info', FUNC_NB)
%    Output: 
%      STR -- function description string of function FUNC_NB
%             FUNC_NB == SSCANF(STR, '%d').  
% Examples:
%    fhs = benchmarks('handles');  
%    f = fh{1}(x);  % evaluates x on the sphere with mod noise, f101
%    f = feval(fh{1}, x); % dito
%    fidx = benchmarks('FunctionIndices');
%
% see also: functions FGENERIC, BENCHMARKS, BENCHMARKSNOISY, BENCHMARKINFOS

% INTERFACE OF BENCHMARK FUNCTIONS
% see benchmarks.m 

%%%-------------------------------------------------------------%%%
function [res, res2] = benchmarksextra(x, strfun, DIM) 
%
  Nfcts = 1;
  baseNum = 200;

  % return valid function IDs (ie numbers)
  if nargin < 1 || ( ...
      ischar(x) && (strcmpi(x, 'FunctionIndices') || strcmpi(x, 'extraFunctionIndices')))
    res = [];
    for i = baseNum + 1:100
      try  % exist does not work under Octave
        eval(['f' num2str(i) '([1; 2]);']);
        res(end+1) = i;
      catch
        if i < baseNum + Nfcts
          disp(sprintf('execution of function %d produces an error', i));
        end
      end
    end 
    if length(res) > Nfcts
      disp(res);
      error([num2str(length(res)) ' > %d functions f1, f2, ..., f24,', ...
             ' f101, ..., f130, f201, ...', ...
             ' are defined. Try to clear the workspace'], Nfcts);
    elseif length(res) < Nfcts
      disp(res);
      error([num2str(length(res)) ' < %d functions are defined. ' ...
             'There might be an execution error '], Nfcts);
    end
  elseif isnumeric(x)
    if isnumeric(strfun)
      res = feval(['f' num2str(strfun)], x);
    else
      res = feval(strfun, x);
    end
  % return function handles
  elseif ischar(x) && strcmpi(x, 'handles')
    res = {};
    for i = 1:100
      try  % exist does not work under Octave
        eval(['f' num2str(i+baseNum) '([1; 2]);']);
        res{i} = str2func(['f' num2str(i+baseNum)]); % eval(['@f' num2str(i)]);
      catch
        if i < Nfcts
          disp(sprintf('execution of function %d produces an error', i+baseNum));
        end
      end
    end
    return;
  elseif ischar(x) && strcmpi(x, 'test')
    error('to be implemented');
    % TODO
    for DIM = [2,3,5,10,20,40]
      % generate test points also outside the domain 
      % pop = generateTestPoints
      for ifun = benchmarksextra('FunctionIndices')
        % evaluate all functions at all the points, using fgeneric
      end
    end
    % dump the result 
  elseif ischar(x) && nargin < 2  % TODO: simplify
    res = benchmarkinfos(); 
    res = res((baseNum+1):end);
    if nargin > 1
      res2 = res;
      res = res2{strfun};
    end
  elseif ischar(x) && strcmpi(x, 'fopt')
    if isnumeric(strfun)
      res = feval(['f' num2str(strfun)], x);
    else
      res = feval(strfun, x);
    end
  elseif ischar(x) 
    res = benchmarkinfos(strfun); 
  else
    error('not a valid case of input');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f201(x, DIM, ntrial)
% sphere with moderate Gauss noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 201; 
  rrseed = 1; 
  aXopt = 0; % actual optimum in input space
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 2;
    end
    x = ones(DIM,1);  % setting all persistent variables
  else
    flginputischar = 0;
  end
  % from here on x is assumed a numeric variable
  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if DIM == 1 
    error('1-D input is not supported');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial; 
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  if isempty(Fopt)
    Fopt = 0;
%     Fopt =1* min(1000, max(-1000, round(100*100*gauss(1,rseed)/gauss(1,rseed+1))/100));
  end 
  Fadd = Fopt;  % value to be added on the "raw" function value
  % DIM-dependent initialization
%   if isempty(lastSize) || lastSize.DIM ~= DIM  
%     Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
%   end
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = x + aXopt;    % shift zero to actual optimum

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2, 1);
  Fval = Ftrue;  % without noise

  %----- FINALIZE -----
  Ftrue = Ftrue + Fadd;
  Fval = Fval + Fadd;

  %----- RETURN INFO ----- 
  if flginputischar  
    if strcmpi(strinput, 'xopt')
      Fval = Fopt;
      Ftrue = Xopt;
    elseif strcmpi(strinput, 'linearTF')
      Fval = Fopt;
      Ftrue = {};
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

% qqq

%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%

function x_opt = computeXopt(seed, DIM)
   % rounded by for digits, but never to zero
   x_opt = 8 * floor(1e4*unif(DIM,seed)')/1e4 - 4;
   idx = x_opt == 0;
   x_opt(idx) = -1e-5;
end

function B = computeRotation(seed, DIM)
% computes an orthogonal basis
  B = reshape(gauss(DIM*DIM,seed), DIM, DIM);
  for i = 1:DIM
    for j = 1:i-1
      B(:,i) = B(:,i) - B(:,i)'*B(:,j) * B(:,j);
    end
    B(:,i) = B(:,i) / sqrt(sum(B(:,i).^2));
  end
end

function g = monotoneTFosc(f)
% maps [-inf,inf] to [-inf,inf] with different constants
% for positive and negative part
   a = 0.1;
   g = f; 
   idx = (f > 0);
   g(idx) = log(f(idx))/a;
   g(idx) = exp(g(idx) + 0.49*(sin(g(idx)) + sin(0.79*g(idx)))).^a;
   idx = (f < 0);
   g(idx) = log(-f(idx))/a;
   g(idx) = -exp(g(idx) + 0.49*(sin(0.55*g(idx)) + sin(0.31*g(idx)))).^a;
end

%---------- pseudo random number generator ------------
function g = gauss(N, seed)
% gauss(N, seed)
% samples N standard normally distributed numbers
% being the same for a given seed
  r = unif(2*N, seed); % in principle we need only half
  g = sqrt(-2*log(r(1:N))) .* cos(2*pi*r(N+1:2*N));
  if any(g == 0)
    g(g==0) = 1e-99;
  end
end

function r = unif(N, inseed)
% unif(N, seed)
%    generates N uniform numbers with starting seed

  % initialization
  inseed = abs(inseed);
  if inseed < 1
    inseed = 1;
  end
  aktseed = inseed;
  for i = 39:-1:0
    tmp = floor(aktseed/127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if aktseed < 0
      aktseed = aktseed + 2147483647;
    end
    if i < 32
      rgrand(i+1) = aktseed;
    end
  end
  aktrand = rgrand(1);

  % sample numbers
  r = zeros(1,N); % makes the function ten times faster(!)
  for i = 1:N
    tmp = floor(aktseed/127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if aktseed < 0
      aktseed = aktseed + 2147483647;
    end
    tmp = floor(aktrand / 67108865);
    aktrand = rgrand(tmp+1);
    rgrand(tmp+1) = aktseed;
    r(i) = aktrand/2.147483647e9;
  end
  if any(r == 0)
    warning('zero sampled(?), set to 1e-99');
    r(r==0) = 1e-99;
  end
end

% for testing and comparing to other implementations, 
%   rand and randn are used only for sampling the noise 
%   should be renamed (ie removed from execution) in the end
%   in 4-D internal rand makes functions 30% faster. 
function res = myrandn(N,M)
   persistent seed
   if isempty(seed)
     seed = 1;
     disp('non-matlab randn');
   else
     seed = seed + 1;  
     if seed > 1e9
       seed = 1;
     end
   end
   res = reshape(gauss(N*M, seed), N, M);
end

function res = myrand(N,M)
   persistent seed
   if isempty(seed)
     seed = 1;
     disp('non-matlab rand');
   else
     seed = seed + 1;
     if seed > 1e9
       seed = 1;
     end
   end
   res = reshape(unif(N*M, seed), N, M); 
end

function Fval = FGauss(Ftrue, beta)
  POPSI = size(Ftrue, 2);
  Fval = Ftrue .* exp(beta * randn(1, POPSI)); % with gauss noise
  TOL = 1e-8; 
  Fval = Fval + 1.01*TOL; 
  idx = Ftrue < TOL;
  Fval(idx) = Ftrue(idx); 
end

function Fval = FUniform(Ftrue, alpha, beta)
  % alpha = 0.49 + 1/DIM;  % alpha * rand must always be smaller than one 
  % beta = 1;              % smaller is easier
  POPSI = size(Ftrue, 2);
  Fval = rand(1,POPSI).^beta .* Ftrue ...
         .* max(1, (10^9 ./ (Ftrue+1e-99)).^(alpha * rand(1, POPSI)));
  TOL = 1e-8; 
  Fval = Fval + 1.01*TOL; 
  idx = Ftrue < TOL;
  Fval(idx) = Ftrue(idx); 
end

function Fval = FCauchy(Ftrue, alpha, p)
  % Cauchy with median 1e3*alpha and with p=0.2, zero otherwise
  % P(Cauchy > 1,10,100,1000) = 0.25, 0.032, 0.0032, 0.00032
  % alpha = 1;
  POPSI = size(Ftrue, 2);
  Fval = Ftrue + alpha * max(0, 1e3 + (rand(1, POPSI) < p) .* ... 
                          randn(1, POPSI) ./ (abs(randn(1, POPSI))+1e-199)); 
  TOL = 1e-8; 
  Fval = Fval + 1.01*TOL; 
  idx = Ftrue < TOL;
  Fval(idx) = Ftrue(idx); 
end

% qqq
%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = template(x, DIM, ntrial)
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = INPUTTHIS; 
  condition = INPUTTHIS;  % for linear transformation
%  alpha = 1;      % 
%  beta = 0.25;       % 
  rrseed = INPUTTHIS;  % original function number
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 2;
    end
    x = ones(DIM,1);  % setting all persistent variables
  else
    flginputischar = 0;
  end
  % from here on x is assumed a numeric variable
  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if DIM == 1 
    error('1-D input is not supported');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial; 
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  if isempty(Fopt)
    Fopt =1* min(1000, max(-1000, round(100*100*gauss(1,rndseed)/gauss(1,rndseed+1))/100));
  end 
  Fadd = Fopt;  % value to be added on the "raw" function value
  % DIM-dependent initialization
  if isempty(lastSize) || lastSize.DIM ~= DIM  
    Xopt =1* computeXopt(rndseed, DIM); % function ID is seed for rotation 
    rotation = computeRotation(rndseed+1e6, DIM); 
    scales = sqrt(condition).^linspace(0, 1, DIM)'; 
    linearTF = diag(scales) * computeRotation(rndseed, DIM); 
    % decouple scaling from function definition
    % linearTF = rotation * linearTF; % or computeRotation(rndseed+1e3, DIM)
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
    % arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI); 
    % arrScales = repmat(scales, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  % x = rotation * x;  % no scaling here, because it would go to the arrExpo
  % x = monotoneTFosc(x); 
  % idx = x > 0;
  % x(idx) = x(idx).^(1 + arrExpo(idx) .* sqrt(x(idx)));  % smooth in zero
  % x = rotation' * x;   % back-rotation
  % x = arrScales .* x;  % scaling, Xopt should be 0 
  x = linearTF * x;    % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = 0; 
  Fval = Ftrue;  % without noise

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 
  Fval = FCauchy(Ftrue, 1, 0.2); 

  %----- FINALIZE -----
  Ftrue = Ftrue + Fadd;
  Fval = Fval + Fadd;

  %----- RETURN INFO ----- 
  if flginputischar  
    if strcmpi(strinput, 'xopt')
      Fval = Fopt;
      Ftrue = Xopt;
    elseif strcmpi(strinput, 'linearTF')
      Fval = Fopt;
      Ftrue = {};
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function