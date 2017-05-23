classdef (Abstract) BbcClientBase < handle
  %BBCCLIENTBASE Summary of this class goes here
  %   Detailed explanation goes here

  properties
  end

  methods (Abstract)
    login(obj);
    numProblems = getNumberOfProblems(obj);
    setProblem(obj, problemID);
    setTrack(obj, trackname);
    dim = getDimension(obj);
    bud = getBudget(obj);
    evals = getEvaluations(obj);
    value = evaluate(obj, point);
    value = safeEvaluate(obj, point, problemID, trackname);
    value = truncate2bounds(obj, value);
  end
end

