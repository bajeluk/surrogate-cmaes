classdef BbcClientShLib < BbcClient
  %BBCCLIENT A client for the Black-Box Competition server.
    
  properties (Access = protected)
    libname
    libhfile
  end % properties

  methods (Access = public)
    function obj = BbcClientShLib(libname, libhfile, username, password, maxTrials)
      % initialization
      obj@BbcClient(username, password, maxTrials)
      
      obj.libname = libname;
      obj.libhfile = libhfile;

      obj.loadBBCompLibrary();
    end

    function cleanup(obj)
      if libisloaded(obj.libname)
        unloadlibrary(obj.libname);
      end
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

    function varargout = call(obj, method, varargin)
      % a wrapper over shared library calls
      nin = length(varargin);

      switch method
        case {'login', 'configure'}
          if nin ~= 2
            obj.throwCallException(2, method, nin);
          else
            val1 = calllib(obj.libname, method, varargin{1}, varargin{2});
            if nargout
              varargout{1} = val1;
            end
          end
        case {'setTrack', 'trackName', 'setProblem'}
          if nin ~= 1
            obj.throwCallException(1, method, nin);
          else
            val1 = calllib(obj.libname, method, varargin{1});
            if nargout
              varargout{1} = val1;
            end
          end
        case {'numberOfTracks', 'numberOfProblems', 'dimension', ...
            'budget', 'evaluations', 'errorMessage'}
          if nin ~= 0
            obj.throwCallException(0, method, nin);
          else
            val1 = calllib(obj.libname, method);
            if nargout
              varargout{1} = val1;
            end
          end
        case 'evaluate'
          if nin ~= 3
            obj.throwCallException(3, method, nin);
          else
            [val1, val2, val3] = calllib(obj.libname, method, ...
              varargin{1}, varargin{2});
            if nargout >= 3
              varargout{3} = val3;
            elseif nargout >= 2
              varargout{2} = val2;
            elseif nargout >= 1
              varargout{1} = val1;
            end
          end
        case 'history'
          if nin ~= 2
            obj.throwCallException(2, method, nin);
          else
            [val1, val2, val3] = calllib(obj.libname, method, ...
              varargin{1}, varargin{2});
            if nargout >= 3
              varargout{3} = val3;
            elseif nargout >= 2
              varargout{2} = val2;
            elseif nargout >= 1
              varargout{1} = val1;
            end
          end
        otherwise
          throw(MException('BbcClient:call', 'Unknown method ''%s''', method));
      end
    end

    function throwCallException(obj, exp_args, method, act_args)
      throw(MException('BbcClient:call', ...
              'Expected %d args in call ''%s'', got %d.', ...
              exp_args, arg_type, method, act_args));
    end
  end % methods

end % classdef

