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
  end

  methods
    function obj = Population(lambda_, dim_)
      % constructor
      obj.lambda = lambda_;
      obj.dim = dim_;
      obj.x = NaN(dim_, lambda_);
      obj.y = NaN(1,    lambda_);
      obj.arx = NaN(dim_, lambda_);
      obj.arz = NaN(dim_, lambda_);
      obj.nPoints = 0;
      obj.origEvaled = false(1, lambda_);
    end

    function obj = addPoints(obj, xNew, yNew, arxNew, arzNew, nOrigEvaled);
      % add new points with their f-values into the population
      % 'nOrigEvaled' - the number of points starting from begiining (index 1)
      %                 which have original fitness value
      nNew = length(yNew);
      if (nNew == 0)
        % nothing to save
        return;
      end
      assert(size(xNew, 2) == nNew, 'Number of points and its y-values are not consistent!');

      obj.x(:, obj.nPoints + [1:nNew])   = xNew;
      obj.y(obj.nPoints + [1:nNew])      = yNew;
      obj.arx(:, obj.nPoints + [1:nNew]) = arxNew;
      obj.arz(:, obj.nPoints + [1:nNew]) = arzNew;
      obj.origEvaled(obj.nPoints + [1:nOrigEvaled]) = true;

      obj.nPoints = obj.nPoints + nNew;
    end

    function [obj, sInd] = sort(obj)
      % sort the points in population according to the f-values
      [sY, sInd] = sort(obj.y);
      obj.y = sY;
      obj.x = obj.x(:, sInd);
      obj.arx = obj.arx(:, sInd);
      obj.arz = obj.arz(:, sInd);
      obj.origEvaled = obj.origEvaled(sInd);
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
