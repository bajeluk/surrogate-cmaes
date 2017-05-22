classdef BbcClient < BbcClientBase
  %BBCCLIENT A client for the Black-Box Competition server.
    
  properties (Access = protected)
    libname
    libhfile
    username
    password
    maxTrials
  end % properties

  methods (Access = public)
    function obj = BbcClient(libname, libhfile, username, password, maxTrials)
      % initialization
      obj.libname = libname;
      obj.libhfile = libhfile;
      obj.username = username;
      obj.password = password;
      obj.maxTrials = maxTrials;

      obj.loadBBCompLibrary();
    end

    function login(obj)
      result = calllib(obj.libname, 'login', obj.username, obj.password);
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

    function safeLogin(obj, delay)
      trial = 1;
      while trial <= obj.maxTrials
        try
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
      result = calllib(obj.libname, 'setTrack', trackname);
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
            obj.safeLogin(maxTrials, 1);
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
      numTracks = calllib(obj.libname, 'numberOfTracks');
      if ~numTracks
        throw(MException('BbcClient:getNumberOfTracks', ...
          'numberOfTracks failed with message: %s', obj.errorMessage()));
      end
    end

    function trackName = getTrackName(obj, i)
      trackName = calllib(obj.libname, 'trackName', i-1); % indexing from 0
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
      numProblems = calllib(obj.libname, 'numberOfProblems');
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
      result = calllib(obj.libname, 'setProblem', problemID);
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
            obj.safeSetTrack(trackname, obj.maxTrials);
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
      dim = calllib(obj.libname, 'dimension');
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
            obj.safeSetProblem(problemID, trackname, obj.maxTrials);
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
      bud = calllib(obj.libname, 'budget');
      if ~bud
        throw(MException('BbcClient:getBudget', ...
          'budget failed with message: %s', obj.errorMessage()));
      end
    end

    function bud = safeGetBudget(obj, problemID, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          bud = obj.getDimension();
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:getBudget')
            warning('getBudget in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetProblem(problemID, trackname, obj.maxTrials);
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
      evals = calllib(obj.libname, 'evaluations');
      if evals < 0
        throw(MException('BbcClient:getEvaluations', ...
          'evaluations failed with message: %s', obj.errorMessage()));
      end
    end

    function evals = safeGetEvaluations(obj, problemID, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          evals = obj.getDimension();
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:getEvaluations')
            warning('getEvaluations in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetProblem(problemID, trackname, obj.maxTrials);
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

    function value = evaluate(obj, point)
      value = libpointer('doublePtr', 1e+100);
      [result, ~, value] = calllib(obj.libname, 'evaluate', point, value);
      if ~result
        throw(MException('BbcClient:evaluate', ...
          'evaluate failed with message: %s', obj.errorMessage()));
      end
    end

    function value = safeEvaluate(obj, point, problemID, trackname)
      trial = 1;
      while trial <= obj.maxTrials
        try
          value = obj.evaluate(point);
          break;
        catch ME
          if strcmp(ME.identifier, 'BbcClient:evaluate')
            warning('evaluate in trial %d / %d failed with message: %s.', ...
              trial, obj.maxTrials, ME.message);
            obj.safeSetProblem(problemID, trackname, obj.maxTrials);
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

    function cleanup(obj)
      if libisloaded(obj.libname)
        unloadlibrary(obj.libname);
      end
    end
  end % methods

  methods (Static)
    function value = truncate2bounds(value)
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
      msg = calllib(obj.libname, 'errorMessage');
      if ~ischar(msg)
        msg = '';
      end
    end
  end % methods
  
end % classdef

