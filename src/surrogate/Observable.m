classdef Observable
  properties
    observers
  end

  methods
    function self = Observable()
      self.observers = {};
    end

    function self = register_observer(self, observer)
      self.observers{end+1} = observer;
    end

    function notify_observers(self, varargin)
      for i = 1:length(self.observers)
        self.observers{i}.notify(self, varargin);
      end
    end
  end
end

