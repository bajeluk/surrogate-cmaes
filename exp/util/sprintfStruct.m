function output = sprintfStruct(s, varargin)
  output = [];
  if length(s) > 1
    for i = 1:length(s)
      output = [output '\\n' sprintfStruct(s(i), 'escape')];
    end
  else
    for fname = fieldnames(s)'
      empty = true;
      if (isnumeric(s.(fname{1})))
        str = num2str(s.(fname{1}));
      elseif (isstr(s.(fname{1})))
        str = s.(fname{1});
      elseif (islogical(s.(fname{1})))
        if (s.(fname{1}))
          str = 'true';
        else
          str = 'false';
        end
      elseif (isstruct(s.(fname{1})) && ~strcmpi(fname{1}, 'parameterSets'))
        str = sprintfStruct(s.(fname{1}), 'escape');
      elseif (iscell(s.(fname{1})))
        c = s.(fname{1});
        str = sprintf('%s ', c{:});
      end
      if (~isempty(str))
        if (nargin > 1 && strcmpi(varargin{1}, 'escape'))
          backN = '\\\\\\n';
        else
          backN = '\\n';
        end
        output = [output sprintf(['%15s: %s' backN], fname{1}, str)];
      end
    end
  end
end

