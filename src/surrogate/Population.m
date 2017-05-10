classdef Population
  % class Population
  %
  % class for carying points in the current population with their
  % fitness values and true/false whether they were evaluated by the
  % original or model fitness
  properties
    lambda
    dim
    x
    y
    arx
    arz
    nPoints
    origEvaled
    isEvaled
    phase
  end

  methods
    function obj = Population(lambda_, dim_)
      % constructor
      obj.lambda = lambda_;
      obj.dim = dim_;
      matrixWidth = 0;
      obj.x = NaN(dim_, matrixWidth);
      obj.y = NaN(1,    matrixWidth);
      obj.arx = NaN(dim_, matrixWidth);
      obj.arz = NaN(dim_, matrixWidth);
      obj.nPoints = 0;
      obj.origEvaled = false(1, matrixWidth);
      obj.phase = (-1) * ones(1, matrixWidth);
      obj.isEvaled = false(1, matrixWidth);
    end

    function obj = addPoints(obj, xNew, yNew, arxNew, arzNew, nOrigEvaled, varargin)
      % add new points with their f-values into the population
      % 'nOrigEvaled' - the number of points starting from begiining (index 1)
      %                 which have original fitness value
      % 'varargin{1}' - the phase number which these points came in
      %                 in DoubleTraineEC:
      %                         Phase 0 -- presample
      %                         Phase 1 -- orig. evaluation of the "best" point(s)
      %                         Phase 2 -- model evaluations of the rest of points
      nNew = size(xNew, 2);
      if (nNew == 0)
        % nothing to save
        return;
      end
      if (isempty(yNew)), yNew = NaN(1, nNew);  end

      % varargin processing
      thisPhase = 0;
      if (nargin >= 7 && ~isempty(varargin{1}))
        thisPhase = varargin{1};
      end
      assert(length(yNew) == nNew, 'Number of points and its y-values are not consistent!');

      % assign the new points' information to the respective fields
      obj.x(:, obj.nPoints + (1:nNew))   = xNew;
      obj.y(obj.nPoints + (1:nNew))      = yNew;
      obj.isEvaled(obj.nPoints + (1:nNew)) = ~isnan(yNew);
      obj.arx(:, obj.nPoints + (1:nNew)) = arxNew;
      obj.arz(:, obj.nPoints + (1:nNew)) = arzNew;
      obj.origEvaled(obj.nPoints + (1:nNew)) = ...
          [true(1,nOrigEvaled) false(1,nNew - nOrigEvaled)];
      obj.phase(obj.nPoints + (1:nNew)) = thisPhase;

      obj.nPoints = obj.nPoints + nNew;
    end

    function obj = updateYValue(obj, xNew, yNew, nOrigEvaled, varargin)
      % Update y-values for the points
      % By default, points which have set   isEvaled == false   are updated,
      % but the points can be specified as vector of indicies into Population
      % in the 'varargin{2}'
      % The 'xNew' parameter is either filled -- in that case the points for which
      % y-values are updated are tested for similarity in X-space, or empty.
      % 'nOrigEvaled' - the number of points starting from yNew(1)
      %                 which have original fitness value
      % 'varargin{1}' - the phase number which these points came in
      %                 in DoubleTraineEC:
      %                         Phase 0 -- presample
      %                         Phase 1 -- orig. evaluation of the "best" point(s)
      %                         Phase 2 -- model evaluations of the rest of points
      % 'varargin{2}' - vector of indices of points into Population which should
      %                 be updated
      nNew = length(yNew);
      if (nNew == 0)
        % nothing to save
        return;
      end

      % varargin processing
      % defalut phase:
      thisPhase = [];
      % default for indices to update: not-evaluated points:
      idxToUpdate = find(~obj.isEvaled(1:obj.nPoints));
      if (nargin >= 6 && ~isempty(varargin{2}))
        idxToUpdate = varargin{2};
      end
      if (nargin >= 5 && ~isempty(varargin{1}))
        thisPhase = varargin{1};
      end
      assert(length(idxToUpdate) >= nNew, 'Number of not-evaluated points and supplied y-values are not consistent!');
      idxToUpdate = idxToUpdate(1:min(nNew, end));

      % check that supplied X's are almost the same as already saved ones
      if (~isempty(xNew))
        isXSimilar = abs(obj.x(:, idxToUpdate) - xNew) < 1000*eps;
        if (~all(all(isXSimilar)))
          warning('The X-values are not similar in %d cases!', sum(any(~isXSimilar,1)));
        end
      end

      % update the y's and other fields
      obj.y(idxToUpdate)      = yNew;
      obj.origEvaled(idxToUpdate(1:nOrigEvaled)) = true;
      if (~isempty(thisPhase))
        obj.phase(idxToUpdate) = thisPhase;
      end
      obj.isEvaled(idxToUpdate) = ~isnan(yNew);
    end

    function [X, Z] = getNotEvaledX(obj)
      % returns matrix of size  lambda x n  with not-so-far evaluated
      % points (ie. which have isEvaled == false)
      X = obj.x(:, ~obj.isEvaled);
      Z = obj.arz(:, ~obj.isEvaled);
    end

    function [X, Z] = getNotOrigEvaledX(obj)
      % returns matrix of size  lambda x n  with points
      % not evaluated by the original fitness (ie. which
      % have isEvaled == false)
      X = obj.x(:, ~obj.origEvaled);
      Z = obj.arz(:, ~obj.origEvaled);
    end

    function X = getOriginalX(obj)
      X = obj.x(:, obj.origEvaled);
    end

    function y = getOriginalY(obj)
      y = obj.y(obj.origEvaled);
    end

    function X = getModelX(obj)
      X = obj.x(:, (obj.isEvaled & ~obj.origEvaled));
    end

    function y = getModelY(obj)
      y = obj.y((obj.isEvaled & ~obj.origEvaled));
    end

    function [obj, xRemoved] = removeNotOrigEvaluated(obj, n, varargin)
      % remove specified number of not original-evaluated points
      % (ie. which have origEvaled == false)
      % if bool vector is given in varargin{1}, the points
      % with true value are removed (should hold n == sum(varargin{1}))
      idxNotOrigEvaled = find(~obj.origEvaled);
      if (nargin >= 3 && ~isempty(varargin{1}))
        idxToRemove = idxNotOrigEvaled(varargin{1});
      else
        idxToRemove = idxNotOrigEvaled;
      end
      if (n > length(idxToRemove))
        warning('removeNotEvaluated(): # of points to be removed is higher than # of not evaluated points. Only available/specified will be deleted.');
      end

      n = min(n, length(idxToRemove));
      idxToRemove = idxToRemove(1:n);

      xRemoved = obj.x(:, idxToRemove);
      obj.x(:, idxToRemove)     = [];
      obj.y(idxToRemove)        = [];
      obj.isEvaled(idxToRemove) = [];
      obj.arx(:, idxToRemove)   = [];
      obj.arz(:, idxToRemove)   = [];
      obj.origEvaled(idxToRemove) = [];
      obj.phase(idxToRemove)    = [];

      obj.nPoints = obj.nPoints - n;
    end

    function [obj, sInd] = sort(obj)
      % sort the points in population according to the f-values
      [sY, sInd] = sort(obj.y);
      obj.y = sY;
      obj.x = obj.x(:, sInd);
      obj.arx = obj.arx(:, sInd);
      obj.arz = obj.arz(:, sInd);
      obj.origEvaled = obj.origEvaled(sInd);
      obj.phase = obj.phase(sInd);
    end

    function fmin = getMinModeled(obj)
      % return the minial non-original fitness
      fmin = min(obj.y(~obj.origEvaled));
    end

    function obj = shiftY(obj, shift)
      % shift the f-values by 'shift'
      obj.y = obj.y + shift;
    end
  end
end
