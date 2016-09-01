classdef FinalPopulation
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
    function obj = FinalPopulation(lambda_, dim_)
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
      [sY, sInd] = sort(obj.y);
      obj.y = sY;
      obj.x = obj.x(:, sInd);
      obj.arx = obj.arx(:, sInd);
      obj.arz = obj.arz(:, sInd);
      obj.origEvaled = obj.origEvaled(sInd);
    end

    function fmin = getMinModeled(obj)
      fmin = min(obj.y(~obj.origEvaled));
    end

    function obj = shiftY(obj, shift)
      obj.y = obj.y + shift;
    end
  end
end
