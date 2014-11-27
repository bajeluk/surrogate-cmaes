classdef GenerationEC
  properties
    origGenerations;
    modelGenerations;
    currentMode         = 'original';
    currentGeneration   = 1;
    lastOriginalGeneration = -1;
    remaining           = 2;
  end

  methods
    function obj = GenerationEC(origGenerations, modelGenerations, varargin)
      % initGenerations' default value 0
      if (nargin >= 3)
        initGenerations = varargin{1};
        currentMode = 'initial';
        remaining = initGenerations;
      else
        initGenerations = 0;
        currentMode = 'original';
        remaining = origGenerations;
      end

      obj.origGenerations = origGenerations;
      obj.modelGenerations = modelGenerations;
      obj.currentGeneration   = 1;
    end

    function result = evaluateOriginal(obj)
      % test whether evalute with the original function
      result = strcmp(currentMode, {'original', 'initial'});
    end

    function result = evaluateModel(obj)
      % test whether evalute with a model
      result = strcmp(currentMode, 'model');
    end 

    function next(obj)
      % change the currentMode if all the generations from
      % the current mode have passed
      obj.remaining = obj.remaining - 1;
      switch obj.currentMode
        case 'initial'
          if (obj.remaining == 0)
            obj.currentMode = 'original';
            obj.remaining = obj.origGenerations;
          end
          obj.lastOriginalGeneration = obj.currentGeneration;
        case 'original'
          if (obj.remaining == 0)
            obj.currentMode = 'model';
            obj.remaining = obj.modelGenerations;
          end
          obj.lastOriginalGeneration = obj.currentGeneration;
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

    function holdOn(obj)
      % call this instead of next() if you want to
      % leave the current mode
      obj.currentGeneration = obj.currentGeneration + 1;
    end

    function gen = getLastOriginalGeneration(obj)
      % get the number of the last generation when the original
      % model was used
      % TODO: add a parameter how many generations to return
      gen = obj.lastOriginalGeneration;
    end
  end
end
