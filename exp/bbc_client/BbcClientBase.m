classdef (Abstract) BbcClientBase
  %BBCCLIENTBASE Summary of this class goes here
  %   Detailed explanation goes here

  properties
  end

  methods (Abstract)
    login(obj);
    numProblems = getNumberOfProblems(obj);
    setProblem(obj, problemID);
    setTrack(obj, trackname);
    getDimension(obj);
    getBudget(obj);
    getEvaluations(obj);
    evaluate(obj, point);
    safeEvaluate(obj, point, problemID, trackname)
  end
end

