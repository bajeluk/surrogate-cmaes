classdef (Abstract) BbcClientBase < handle
  %BBCCLIENTBASE Summary of this class goes here
  %   Detailed explanation goes here

  properties
    dim
  end

  methods (Abstract)
    configure(obj, history, logfilepath);
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
    [point, value] = history(obj, index);
  end
end

