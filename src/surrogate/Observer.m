classdef Observer
%OBSERVER -- Observer design pattern
%
% Decoupled data/model and view component. Observer is an independent
% object which can be registred by (multiple) subjects (Observables' subclass)
% in order to view/log information from subjects. These subjects send 
% either specific kind of information (in varargin) to be viewed, or the whole
% subject itself and Observer finds the right information what it needs
%
% subject = Observable();
% observer = Observer();
% subject = observer.registerObservable(subject)
% subject.notify_observers('test')
  methods
    function self = Observer()
    end

    function observable = registerObservable(self, observable)
      observable = observable.register_observer(self);
    end

    function notify(self, observable, varargin)
      disp('Got:');
      if (nargin > 2)
        disp(varargin{1});
      end
    end
  end
end
