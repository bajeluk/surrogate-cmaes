classdef BbcClientTcp < BbcClient
  %BBCCLIENT A client for the Black-Box Competition server.

  properties (Access = protected)
    hostname
    port
    timeout
    connectTimeout
    tcpclient
  end % properties

  methods (Access = public)
    function obj = BbcClientTcp(hostname, port, username, password, ...
        timeout, connectTimeout, maxTrials, delay)
      % initialization
      obj@BbcClient(username, password, maxTrials, delay);

      obj.hostname = hostname;
      obj.port = port;
      obj.username = username;
      obj.password = password;
      obj.maxTrials = maxTrials;
      obj.timeout = timeout;
      obj.connectTimeout = connectTimeout;

      try
        obj.tcpclient = tcpclient(obj.hostname, obj.port, ...
          'Timeout', timeout, 'ConnectTimeout', connectTimeout);
      catch ME
        throw(MException('BbcClient:connect', ...
          'Could not connect to proxy server %s:%d.', ...
          obj.hostname, obj.port));
      end
    end

    function cleanup(obj)
      % tcp connection will be closed upon object removal
    end
  end % methods

  methods (Access = protected)
    function varargout = call(obj, method, varargin)
      % a wrapper over TCP calls
      nin = 0;
      ignored_types = {'lib.pointer'};
      
      for i = 1:length(varargin)
        if ~ismember(class(varargin{i}), ignored_types)
          nin = nin + 1;
        end
      end

      switch method
        case {'login', 'configure'}
          if nin ~= 2
            obj.throwCallException(2, method, nin);
          end
        case {'setTrack', 'trackName', 'setProblem', 'evaluate', 'history'}
          if nin ~= 1
            obj.throwCallException(1, method, nin);
          end
        case {'numberOfTracks', 'numberOfProblems', 'dimension', ...
          'budget', 'evaluations', 'errorMessage'}
          if nin ~= 0
            obj.throwCallException(0, method, nin);
          end
        otherwise
          throw(MException('BbcClient:call', 'Unknown method ''%s''', method));
      end

      % format input arguments into string
      strargin = cell(1, nin);

      for i = 1:nin
        arg = varargin{i};
        switch class(arg)
          case 'char'
            if ~isempty(strfind(arg, '"'))
              throw(MException('BbcClient:call', ...
                'Quotes inside string arguments not allowed.'));
            end
            strargin{i} = ['"', arg, '"'];
          case {'int8', 'int16', 'int32', 'int64', ...
              'uint8', 'uint16', 'uint32', 'uint64'}
            strarg = arrayfun(@(x) num2str(x), arg, ...
              'UniformOutput', false);
            if length(strarg) == 1
              strargin{i} = strarg{:};
            else
              strargin{i} = ['[', strjoin(strarg, ' '), ']'];
            end
          case 'double'
            strarg = arrayfun(@(x) num2str(x, '%e'), arg, ...
              'UniformOutput', false);
            if length(strarg) == 1
              strargin{i} = strarg{:};
            else
              strargin{i} = ['[', strjoin(strarg, ' '), ']'];
            end
          case ignored_types
            % ignore
          otherwise
            throw(MException('BbcClient:call', ...
              'Unsupported type ''%s''', class(arg)));
        end
      end

      request = sprintf('%s(%s)', method, strjoin(strargin, ','));

      try
        write(obj.tcpclient, [uint8(request) uint8(10)]);
        
        response = [];
        while isempty(response) || response(end) ~= 10
          data = read(obj.tcpclient);
          response = [response data];
        end

        response = strsplit(char(response(1:end-1)), ',');

        if strcmpi(response{1}, 'ERROR')
          throw(MException('BbcClient:call', ...
            'Protocol error: %s', response{2}));
        end
      catch ME
        throw(MException('BbcClient:call', ...
          'TCP request failed with error %s, message: %s', ME.identifier, ME.message));
      end

      switch method
        case {'login', 'configure', 'setTrack', 'setProblem', 'numberOfTracks', ...
          'numberOfProblems', 'dimension', 'budget', 'evaluations'}
          if nargout
            varargout{1} = str2double(response{1});
          end
        case {'errorMessage', 'trackName'}
          if nargout
            varargout{1} = strrep(response{1}, '"', '');
          end
        case 'evaluate'
          if nargout >= 3
            if length(response) < 2
              varargout{3} = NaN;
            else
              varargout{3} = str2double(response{2});
            end
          end

          if nargout >= 2
            varargout{2} = [];
          end

          if nargout >= 1
            varargout{1} = str2double(response{1});
          end
        case 'history'
          varargout{1} = str2double(response{1});
          varargout{2} = str2num(response{2});
          varargout{3} = str2double(response{3});
      end
    end

    function throwCallException(obj, exp_args, method, act_args)
      throw(MException('BbcClient:call', ...
              'Expected %d args in call ''%s'', got %d.', ...
              exp_args, arg_type, method, act_args));
    end
  end % methods

end % classdef

