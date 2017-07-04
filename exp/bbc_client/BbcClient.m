classdef BbcClient < BbcClientBase
  %BBCCLIENT A client for the Black-Box Competition server.

  properties (Access = protected)
    username
    password
    maxTrials
  end % properties

  methods (Access = public)
    function obj = BbcClient(username, password, maxTrials)
      % initialization
      obj.username = username;
      obj.password = password;
      obj.maxTrials = maxTrials;
      obj.dim = NaN;
    end

    function login(obj)
      result = obj.call('login', obj.username, obj.password);
      if ~result
        msg = obj.errorMessage();
        if strcmp(msg, 'already logged in')
          return
        else
          throw(MException('BbcClient:login', ...
            'login failed with message: %s', msg));
        end
      end
    end

    function configure(obj, history, logfilepath)
      result = obj.call('configure', uint8(history), logfilepath);
      if ~result
        throw(MException('BbcClient:configure', ...
          'configure failed with message: %s', obj.errorMessage()));
      end
    end

    function [point, value] = history(obj, index)
      value = libpointer('doublePtr', 1e+100);
      point = libpointer('doublePtr', zeros(1,obj.dim));
      [result, point, value] = obj.call('history', uint16(index), point, value);

      if ~result
        throw(MException('BbcClient:history', ...
          'history failed with message: %s', obj.errorMessage()));
      end
    end

    function safeLogin(obj, delay)
      trial = 1;
      while trial <= obj.maxTrials
        try
%           fprintf('Debug: login trial %d / %d\n', ...
%             trial, obj.maxTrials);
          obj.login();
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:login')
            warning('Login in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            pause(delay);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeLogin', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function setTrack(obj, trackname)
      result = obj.call('setTrack', trackname);
      if ~result
        throw(MException('BbcClient:setTrack', ...
          'setTrack failed with message: %s', obj.errorMessage()));
      end
    end
    
    function safeSetTrack(obj, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          obj.setTrack(trackname);
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:setTrack')
            warning('setTrack in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeLogin(1);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeSetTrack', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function numTracks = getNumberOfTracks(obj)
      numTracks = obj.call('numberOfTracks');
      if ~numTracks
        throw(MException('BbcClient:getNumberOfTracks', ...
          'numberOfTracks failed with message: %s', obj.errorMessage()));
      end
    end

    function trackName = getTrackName(obj, i)
      trackName = obj.call('trackName', int8(i-1)); % indexing from 0
      if isempty(trackName)
        throw(MException('BbcClient:getTrackName', ...
          'trackName failed with message: %s', obj.errorMessage()));
      end
    end

    function trackNames = getTrackNames(obj)
      try
        numTracks = obj.getNumberOfTracks();
        trackNames = cell(1, numTracks);

        for i = 1:numTracks
          trackNames{i} = obj.getTrackName(i);
        end
      catch ME
        rethrow(ME);
      end
    end

    function numProblems = getNumberOfProblems(obj)
      numProblems = obj.call('numberOfProblems');
      if ~numProblems
        throw(MException('BbcClient:getNumberOfProblems', ...
          'numberOfProblems failed with message: %s', obj.errorMessage()));
      end
    end

    function numProblems = safeGetNumberOfProblems(obj)
      trial = 1;
      while trial <= obj.maxTrials
        try
          numProblems = obj.getNumberOfProblems();
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:getNumberOfProblems')
            warning('getNumberOfProblems in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeLogin(obj.maxTrials, 1);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeGetNumberOfProblems', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function setProblem(obj, problemID)
      result = obj.call('setProblem', uint16(problemID-1)); % indexing from 0
      if ~result
        throw(MException('BbcClient:setProblem', ...
          'setProblem failed with message: %s', obj.errorMessage()));
      end
    end

    function safeSetProblem(obj, problemID, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          obj.setProblem(problemID);
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:setProblem')
            warning('setProblem in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetTrack(trackname);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeSetProblem', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function dim = getDimension(obj)
      dim = obj.call('dimension');
      obj.dim = dim;
      if ~dim
        throw(MException('BbcClient:getDimension', ...
          'dimension failed with message: %s', obj.errorMessage()));
      end
    end

    function dim = safeGetDimension(obj, problemID, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          dim = obj.getDimension();
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:getDimension')
            warning('getDimension in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetProblem(problemID, trackname);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeGetDimension', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function bud = getBudget(obj)
      bud = obj.call('budget');
      if ~bud
        throw(MException('BbcClient:getBudget', ...
          'budget failed with message: %s', obj.errorMessage()));
      end
    end

    function bud = safeGetBudget(obj, problemID, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          bud = obj.getBudget();
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:getBudget')
            warning('getBudget in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetProblem(problemID, trackname);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeGetBudget', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function evals = getEvaluations(obj)
      evals = obj.call('evaluations');
      if evals < 0
        throw(MException('BbcClient:getEvaluations', ...
          'evaluations failed with message: %s', obj.errorMessage()));
      end
    end

    function evals = safeGetEvaluations(obj, problemID, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          evals = obj.getEvaluations();
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:getEvaluations')
            warning('getEvaluations in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetProblem(problemID, trackname);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeGetEvaluations', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function value = evaluate(obj, point, varargin)
      if ~isempty(varargin) && varargin{1}
        point = obj.truncate2bounds(point);
      end

      value = libpointer('doublePtr', 1e+100);
      [result, ~, value] = obj.call('evaluate', point, value);

      if ~result
        msg = obj.errorMessage();
        if strcmp(msg, 'evaluation budget exceeded')
          value = Inf;
        else
          throw(MException('BbcClient:getEvaluations', ...
            'evaluations failed with message: %s', obj.errorMessage()));
        end
      end
    end

    function value = safeEvaluate(obj, point, problemID, trackname, varargin)
      trial = 1;
      while trial <= obj.maxTrials
        try
          value = obj.evaluate(point, varargin{:});
          break;
        catch ME
          fields = strsplit(ME.identifier, ':');
          if strcmp(fields{1}, 'BbcClient')
            warning('evaluate in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetProblem(problemID, trackname);
            trial = trial + 1;
          else
            rethrow(ME);
          end
        end
      end

      if trial > obj.maxTrials
        throw(MException('BbcClient:safeGetEvaluations', ...
          'Exhausted maximum number of trials: %d', obj.maxTrials));
      end
    end

    function value = truncate2bounds(~, value)
      value = max(0.0, min(1.0, value));
    end
  end % methods

  methods (Access = protected)
    function loadBBCompLibrary(obj)
      if libisloaded(obj.libname) % first unload if already loaded
        unloadlibrary(obj.libname);
      end

      loadlibrary(obj.libname, obj.libhfile);

      libfuns = libfunctions(obj.libname, '-full');
      if ~numel(libfuns)
        throw(MException('BbcClient:loadlibrary', ['Error loading library.' ...
        ' Is library in path?']));
      end
    end

    function msg = errorMessage(obj)
      msg = obj.call('errorMessage');
      if ~ischar(msg)
        msg = '';
      end
    end
  end % methods

  methods (Abstract, Access = public)
    cleanup(obj);
  end
  
  methods (Abstract, Access = protected)
    varargout = call(obj, method, varargin);
  end % methods

end % classdef

