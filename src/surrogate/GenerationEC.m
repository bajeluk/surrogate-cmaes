classdef GenerationEC < handle
  properties
    origGenerations;
    modelGenerations;
    currentMode         = 'original';
    currentGeneration   = 1;
    lastOriginalGenerations = [];
    remaining           = 2;
  end

  methods
    function obj = GenerationEC(origGenerations, modelGenerations, varargin)
      % initGenerations' default value 0
      if (nargin >= 3  &&  varargin{1} > 0)
        initGenerations = varargin{1};
        obj.currentMode = 'initial';
        obj.remaining = initGenerations;
      else
        initGenerations = 0;
        obj.currentMode = 'original';
        obj.remaining = origGenerations;
      end

      obj.origGenerations = origGenerations;
      obj.modelGenerations = modelGenerations;
      obj.currentGeneration   = 1;
    end

    function result = evaluateOriginal(obj)
      % test whether evalute with the original function
      result = any(strcmp(obj.currentMode, {'original', 'initial'}));
    end

    function result = isNextOriginal(obj)
      % check whether there will be 'original' mode after calling next()
      result = (any(strcmp(obj.currentMode, {'original', 'initial'})) && (obj.remaining > 1)) ...
          || (strcmp(obj.currentMode, 'model') && obj.remaining == 1);
    end

    function result = evaluateModel(obj)
      % test whether evalute with a model
      result = strcmp(obj.currentMode, 'model');
    end 

    function obj = next(obj)
      % change the currentMode if all the generations from
      % the current mode have passed
      obj.remaining = obj.remaining - 1;
      switch obj.currentMode
        case 'initial'
          if (obj.remaining == 0)
            obj.currentMode = 'original';
            obj.remaining = obj.origGenerations;
          end
          obj.lastOriginalGenerations = [obj.lastOriginalGenerations obj.currentGeneration];
        case 'original'
          if (obj.remaining == 0)
            obj.currentMode = 'model';
            obj.remaining = obj.modelGenerations;
          end
          obj.lastOriginalGenerations = [obj.lastOriginalGenerations obj.currentGeneration];
        case 'model'
          if (obj.remaining == 0)
            obj.currentMode = 'original';
            obj.remaining = obj.origGenerations;
          end
        otherwise
          error('GenerationEC: wrong currentMode.');
      end
      obj.currentGeneration = obj.currentGeneration + 1;
    end

    function obj = holdOn(obj)
      % call this instead of next() if you want to
      % leave the current mode
      if (any(strcmp(obj.currentMode, {'original', 'initial'})))
        obj.lastOriginalGenerations = [obj.lastOriginalGenerations obj.currentGeneration];
      end
      obj.currentGeneration = obj.currentGeneration + 1;
    end

    function obj = setNextOriginal(obj)
      % set the next generation and currentMode to 'original'
      % later in the same generation, next() is expected to be called
      obj.currentMode = 'original';
      obj.remaining = 2;
    end

    function gens = getLastOriginalGenerations(obj, n)
      % get the numbers of the last n generations when the original
      % model was used
      startID = length(obj.lastOriginalGenerations) - n + 1;
      if (startID <= 0)
        warning('GenerationEC.getLastOriginalGenerations(): not enough data in the archive');
        startID = 1;
      end
      gens = obj.lastOriginalGenerations(startID:end);
    end
  end
end
