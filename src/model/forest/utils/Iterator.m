classdef (Abstract) Iterator < handle
  methods (Abstract)
    r = hasNext(obj)
    % whether next element is available
    r = next(obj)
    % retrieves the next element
  end
end