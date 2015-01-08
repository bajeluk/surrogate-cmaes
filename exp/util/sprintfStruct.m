function output = sprintfStruct(s, varargin)
  output = [];
  for fname = fieldnames(s)'
    empty = true;
    if (isnumeric(s.(fname{1})))
      str = num2str(s.(fname{1}));
    elseif (isstr(s.(fname{1})))
      str = s.(fname{1});
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

