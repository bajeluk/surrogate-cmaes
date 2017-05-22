classdef BbcClientBBOB
  %BBCCLIENTBBOB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods (Access = public)
    function obj = BbcClientBBOB()
    end

    function login(obj)
    end

    function setTrack(obj, trackname)
    end

    function setProblem(obj, problemID)
    end

    function getDimension(obj)
    end

    function getBudget(obj)
    end

    function getEvaluations(obj)
    end

    function evaluate(obj, point)
    end

    function safeEvaluate(obj, point, problemID, trackname)
    end
  end % methods
end

