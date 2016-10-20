classdef OrigRatioUpdaterRankDiff2 < OrigRatioUpdaterAbstractError
% OrigRatioUpdaterRankDiff2 -- concrete subclass with RankDiff error, version 2
%
% Error is calculating based on rankings of model and original f-values using
% the function errRankMu().
%
% 2016-10-19: version 2 corresponds to the modified OrigRatioUpdater which
%             already implements update() method and leaves only computeErr()
%             for sublcasses
  properties
  end

  methods
    function err = computeErr(obj, modelY, origY, varargin)
      err = errRankMu(modelY, origY, varargin{1});
    end

    function obj = OrigRatioUpdaterRankDiff2(ec, parameters)
      % just call predecessor's constructor
      obj = obj@OrigRatioUpdaterAbstractError(ec, parameters);
    end
  end
end
