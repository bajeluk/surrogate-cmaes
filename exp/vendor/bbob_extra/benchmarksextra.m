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
  Nfcts = 2;
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
    res{1} = @f201;
    res{2} = @f202;
    res{3} = @f203;
    res{4} = @f204;
    res{5} = @f205;
    res{6} = @f206;
    res{7} = @f207;
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
% Sum of Squares Clustering Problem
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrScales arrExpo rseed
  persistent dataset minVal maxVal

  funcID = 201;
  rrseed = 1; 
  aXopt = 0; % actual optimum in input space
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 8;
    end
    x = ones(DIM,1);  % setting all persistent variables
  else
    flginputischar = 0;
  end
  % from here on x is assumed a numeric variable
  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if mod (DIM, 4) ~= 0 
    error('Dimension must be divisable by 4');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial; 
    
    dataset = importdata('iris.txt');
    
    rng(ntrial);
    theta = randi(360);
    R = eye(4, 4);
    R(1, 1) = cosd(theta);
    R(1, 2) = sind(theta);
    R(2, 1) = -sind(theta);
    R(2, 2) = cosd(theta);
    dataset = (R*dataset')';
    
    theta = randi(360);
    R = eye(4, 4);
    R(2, 2) = cosd(theta);
    R(2, 3) = sind(theta);
    R(3, 2) = -sind(theta);
    R(3, 3) = cosd(theta);
    dataset = (R*dataset')';
    
    theta = randi(360);
    R = eye(4, 4);
    R(1, 1) = cosd(theta);
    R(1, 3) = -sind(theta);
    R(3, 1) = sind(theta);
    R(3, 3) = cosd(theta);
    dataset = (R*dataset')';
    
    theta = randi(360);
    R = eye(4, 4);
    R(1, 1) = cosd(theta);
    R(1, 4) = sind(theta);
    R(4, 1) = -sind(theta);
    R(4, 4) = cosd(theta);
    dataset = (R*dataset')';
    
    theta = randi(360);
    R = eye(4, 4);
    R(2, 2) = cosd(theta);
    R(2, 4) = -sind(theta);
    R(4, 2) = sind(theta);
    R(4, 4) = cosd(theta);
    dataset = (R*dataset')';
    
    theta = randi(360);
    R = eye(4, 4);
    R(3, 3) = cosd(theta);
    R(3, 4) = -sind(theta);
    R(4, 3) = sind(theta);
    R(4, 4) = cosd(theta);
    dataset = (R*dataset')';
    
    minVal = min(min(dataset));
    maxVal = max(max(dataset));
    switch DIM
        case 8
            Fopt = 152.347951760358;
        case 12
            Fopt = 78.8514414261460;
        case 16
            Fopt = 57.2284732142857;
        case 20
            Fopt = 46.4461820512821;
        case 24
            Fopt = 39.0399872460872;
        case 28
            Fopt = 34.2982296650718;
        case 32
            Fopt = 29.9889439507861;
        case 36
            Fopt = 27.7860924173081;
        case 40
            Fopt = 25.834054819972511;
        otherwise
            Fopt = 0;
            disp('Fopt is unknown. Using Fopt = 0');
    end
    
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
  end 
  
  x = normalize_arr(x, minVal, maxVal);

  
  res = [];
  for i = 1:POPSI
      res = [res fitnessclustsse(x(:, i), dataset)];
  end
  
  %----- COMPUTATION core -----
  Ftrue = res;
  Fval = res;  % without noise

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

function [Fval, Ftrue] = f202(x, DIM, ntrial)
% Sum of Squares Clustering Problem
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrScales arrExpo rseed
  persistent dataset minVal maxVal

  funcID = 202;
  rrseed = 1; 
  aXopt = 0; % actual optimum in input space
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 4;
    end
    x = ones(DIM,1);  % setting all persistent variables
  else
    flginputischar = 0;
  end
  % from here on x is assumed a numeric variable

  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if mod (DIM, 2) ~= 0 
    error('Dimension must be divisable by 2');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial; 
    
    dataset = importdata('ruspini.txt');
    
    %Rotate 
    rng(ntrial);
    theta = randi(360);
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    points = dataset';
    dataset = (R*points)';
    minVal = min(min(dataset));
    maxVal = max(max(dataset));
    
    %Known bests
    switch DIM
        case 4
            Fopt = 89337.832142857;
        case 6
            Fopt = 51063.4750456705;
        case 8
            Fopt = 12881.0512361466;
        case 10
            Fopt = 10126.7197881828;
        case 12
            Fopt = 8575.40687645688;
        case 14
            Fopt = 7126.19854312354;
        case 16
            Fopt = 6149.639019314019;
        case 18
            Fopt = 5181.65183982684;
        case 20
            Fopt = 4446.28214285714;
        otherwise
            Fopt = 0;
            disp('Fopt is unknown. Using Fopt = 0');
    end
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  if isempty(Fopt)
    Fopt = 0;
  end 
  
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
  end
  x = normalize_arr(x, minVal, maxVal);
  
  res = [];
  for i = 1:POPSI
      res = [res fitnessclustsse(x(:, i), dataset)];
  end
  
  %----- COMPUTATION core -----
  Ftrue = res;
  Fval = res;  % without noise

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

end

function [Fval, Ftrue] = f203(x, DIM, ntrial)
% Sum of Squares Clustering Problem
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrScales arrExpo rseed
  persistent dataset minVal maxVal

  funcID = 203;
  rrseed = 1; 
  aXopt = 0; % actual optimum in input space
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 6;
    end
    x = ones(DIM,1);  % setting all persistent variables
  else
    flginputischar = 0;
  end
  % from here on x is assumed a numeric variable
  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if mod (DIM, 3) ~= 0 
    error('Dimension must be divisable by 3');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial; 
    
    dataset = importdata('german_postal.txt');
    
    rng(ntrial);
    theta = randi(360);
    
    R = eye(3, 3);
    R(2, 2) = cosd(theta);
    R(2, 3) = -sind(theta);
    R(3, 2) = sind(theta);
    R(3, 3) = cosd(theta);
    dataset = (R*dataset')';
    
    R = eye(3, 3);
    R(1, 1) = cosd(theta);
    R(3, 1) = -sind(theta);
    R(1, 3) = sind(theta);
    R(3, 3) = cosd(theta);
    dataset = (R*dataset')';
    
    R = eye(3, 3);
    R(1, 1) = cosd(theta);
    R(1, 2) = -sind(theta);
    R(2, 1) = sind(theta);
    R(2, 2) = cosd(theta);
    dataset = (R*dataset')';
      
    
    minVal = min(min(dataset));
    maxVal = max(max(dataset));
    
    switch DIM
        case 6
            Fopt = 6.025472220938822e+11;
        case 9
            Fopt = 2.94506562778027e+11;
        case 12
            Fopt = 1.04474664100716e+11;
        case 15
            Fopt = 5.97615267205269e+10;
        case 18
            Fopt = 3.59085384380298e+10;
        case 21
            Fopt = 2.19832076153985e+10;
        case 24
            Fopt = 1.33854150525813e+10;
        case 27
            Fopt = 8.423750573163314e+09;
        case 30
            Fopt = 6.44648364287013e+09;
        otherwise
            Fopt = 0;
            disp('Fopt is unknown. Using Fopt = 0');
    end
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  if isempty(Fopt)
    Fopt = 0;
  end 
  
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
  end 
  
  x = normalize_arr(x, minVal, maxVal);
  
  res = [];
  for i = 1:POPSI
      res = [res fitnessclustsse(x(:, i), dataset)];
  end
  
  %----- COMPUTATION core -----
  Ftrue = res;
  Fval = res;  % without noise

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

end

function [Fval, Ftrue] = f204(x, DIM, ntrial)
% Neural network learning Problem
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrScales arrExpo rseed
  persistent dataset

  funcID = 204;
  rrseed = 1; 
  aXopt = 0; % actual optimum in input space
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 10;
    end
    x = ones(DIM,1);  % setting all persistent variables
  else
    flginputischar = 0;
  end
  % from here on x is assumed a numeric variable
  x = normalize_arr(x, -100, 100);
  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if mod ((DIM - 1), 3) ~= 0 
    error('Dimension must be according to 3k + 1');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial; 
    dataset = f204Data(ntrial);
    bestResult = 9000;
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  if isempty(Fopt)
    Fopt = 0;
  end 
  
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
  end 
  
  res = zeros(1, POPSI);
  for i = 1:POPSI
      params = x(:, i);
      for j = 1:length(dataset)
        data = dataset(j, :);
        prediction = params(1);
        for k = 0:(((DIM - 1) / 3) - 1)
            prediction = prediction + params(k*3+2) * tanh(params(k*3+3) + (data(1) * params(k*3+4)));
        end
        res(i) = res(i) + ((prediction - data(2)) ^ 2);
      end
  end
  res = res / length(dataset);
  
  %----- COMPUTATION core -----
  Ftrue = res;
  Fval = res;  % without noise

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
end

function [Fval, Ftrue] = f205(x, DIM, ntrial)
% Sammon mapping Problem
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrScales arrExpo rseed
  persistent dataset orig_dist data

  funcID = 205;
  rrseed = 1; 
  aXopt = 0; % actual optimum in input space
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 4;
    end
    x = ones(DIM,1);  % setting all persistent variables
  else
    flginputischar = 0;
  end
  % from here on x is assumed a numeric variable
  x = normalize_arr(x, -50, 50);
  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if mod (DIM, 2) ~= 0 
    error('Dimension must be divisible by 2');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial;
    if (ntrial < 250)
        data = importdata('virus3.txt');
    else
        data = importdata('samp05.txt');
    end
    s = rng;
    rng(ntrial);
    dataset = datasample(data, DIM / 2, 'Replace', false);  
    rng(s);
    orig_dist = pdist(dataset);
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  if isempty(Fopt)
    Fopt = 0;
  end 
  
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM;
  end 
  
  res = zeros(1, POPSI);
  
  for i = 1:POPSI
    index = 1;
    solution = x(:, i);
    sol_dist = pdist(vec2mat(solution, 2));
    for j = 1:(DIM / 2)
        for k = (j+1):(DIM / 2)
            res(i) = res(i) + ((sol_dist(index) - orig_dist(index)).^2) / (orig_dist(index).^2);
            index = index + 1;
        end
    end
  end
 
  %----- COMPUTATION core -----
  Ftrue = res;
  Fval = res;  % without noise

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
end

function [Fval, Ftrue] = f206(x, DIM, ntrial)
% Coordinates on map Problem
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrScales arrExpo rseed
  persistent dataset

  funcID = 206;
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
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    rseed = rrseed + 1e4 * ntrial; 
    if ntrial == 1
        dataset = importdata('earthelevation.txt');
    else
        data = load('detailed_evelation.mat');
        dataset = imresize(data.data, 1 / (ntrial - 1));
    end
    Fopt = max(max(dataset));    
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  
  % from here on x is assumed a numeric variable
  [sizeX, sizeY] = size(dataset);
  x(1, :) = normalize_arr(x(1, :), 1, sizeY);
  x(2, :) = normalize_arr(x(2, :), 1, sizeX);
  [DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
  if DIM ~= 2 
    error('Dimension must be 2');
  end
  
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
  end
  
  x = round(x);
  res = zeros(1, POPSI);
  for i = 1:POPSI
      coordinates = x(:, i);
      res(i) = Fopt - dataset(coordinates(2), coordinates(1));
  end
  
  %----- COMPUTATION core -----
  Ftrue = res;
  Fval = res;  % without noise

  %----- RETURN INFO ----- 
  if flginputischar  
    if strcmpi(strinput, 'xopt')
      Fval = 0;
      Ftrue = Xopt;
    elseif strcmpi(strinput, 'linearTF')
      Fval = 0;
      Ftrue = {};
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = 0;
    end
  end
end

function [Fval, Ftrue] = f207(x, DIM, ntrial)
% Buoy Placement Problem
% ntrial determines configuration of simulator
% ntrial format is following: frwh
% f = frequencies = 1 - 2 - 3 - 5 - 10 - 25
% r = radius = 2.0 - 2.5 - 3.2 - 4.0 - 5.0
% w = width = 10 - 25 - 50 - 100 - 200 - 350 - 500 - 1000 - 2500
% h = height = 10 - 25 - 50 - 100 - 200 - 350 - 500 - 1000 - 2500
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrScales arrExpo rseed
  persistent dataset
  persistent frequencies radius width height

  funcID = 207;
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
  if mod (DIM, 2) ~= 0 
    error('Dimension must be divisable by 2');
  end
  
  %----- INITIALIZATION -----
  if nargin > 2     % set seed depending on trial index
    Fopt = [];      % clear previous settings for Fopt
    lastSize = [];  % clear other previous settings
    tmp = ntrial;
    w_index = mod(tmp, 10);
    tmp = tmp - w_index;
    h_index = mod(tmp, 100) / 10;
    tmp = tmp - h_index * 10;
    r_index = mod(tmp, 1000) / 100;
    tmp = tmp - r_index * 100;
    f_index = mod(tmp, 10000) / 1000;
    tmp = [1 2 3 5 10 25];
    frequencies = tmp(f_index);
    tmp = [2.0 2.5 3.2 4.0 5.0];
    radius = tmp(r_index);
    tmp = [10 25 50 100 200 350 500 1000 2500];
    width = tmp(w_index);
    tmp = [10 25 50 100 200 350 500 1000 2500];
    height = tmp(h_index);
  elseif isempty(rseed)
    rseed = rrseed; 
  end
  if isempty(Fopt)
    Fopt = 1;
  end 
  
  Xopt = zeros(DIM, 1);
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
  end
  
  x(1, :) = normalize_arr(x(1, :), 1, width);
  x(2, :) = normalize_arr(x(2, :), 1, height);
  x = round(x);
  res = zeros(1, POPSI);
  buoys_count = DIM / 2;
  for i = 1:POPSI
      spheres_x = x(1:buoys_count, i);
      spheres_y = x(buoys_count+1:DIM, i);
      radiusSphere = radius*ones(1, buoys_count);
      
      try
        if (frequencies == 1)
            [MlRCWWeighted, MlParrayAverage] = arrayBuoyPlacementSingleFreq(radiusSphere, spheres_x, spheres_y, 0.7, 0);
        else
            [MlRCWWeighted, MlParrayAverage] = arrayBuoyPlacement_v2(radiusSphere, spheres_x, spheres_y, frequencies);
            end
        res(i) = Fopt - MlRCWWeighted;
      catch ME
        %Simulation can throw exception -> 2 buoyes are on the same
        %location
        res(i) = Fopt;
      end
  end
  
  %----- COMPUTATION core -----
  Ftrue = res;
  Fval = res;  % without noise

  %----- RETURN INFO ----- 
  if flginputischar  
    if strcmpi(strinput, 'xopt')
      Fval = 0;
      Ftrue = Xopt;
    elseif strcmpi(strinput, 'linearTF')
      Fval = 0;
      Ftrue = {};
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = 0;
    end
  end
end

function normalized = normalize_arr(array, x, y)
     % Normalize to [0, 1]:
     m = -5;
     range = 5 - m;
     array = (array - m) / range;

     % Then scale to [x,y]:
     range2 = y - x;
     normalized = (array*range2) + x;
end

function res = f204Data(ntrial)
    x = importdata('x.txt');
    switch ntrial
        case 1
            res = [x x.^2];
        case 2
            res = [x sin(x)];
        case 3
            res = [x abs(x)];
        case 4
            res = [x sinh(x)];
        case 5
            res = [x cos(x)];
        case 6
            res = [x tan(x)];
        case 7
            res = [x cosh(x)];
        case 8
            res = [x tanh(x)];
        case 9
            res = [x cot(x)];
        case 10
            res = [x coth(x)];
    end
            
end

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