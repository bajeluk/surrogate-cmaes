classdef CMCell
  properties 
    center % cell center
    dim    % data dimension
    lb     % lower bound
    ub     % upper bound
    minx   % point of minimal value
    maxx   % point of maximal value
    miny   % minimal y
    maxy   % maximal y
    X      % data in cell
    y      % data values in cell
  end
  
  methods
    function obj = CMCell(X, y, settings)
      % constructor
      if nargin < 3
        if nargin < 1
          X = [];
          y = [];
        end
        settings = struct();
      end
      obj.X = X;
      obj.y = y;
      obj.dim = defopts(settings, 'dim', size(X, 2));
      obj.lb = defopts(settings, 'lb', min(X));
      obj.ub = defopts(settings, 'ub', max(X));
      obj.center = obj.lb + (obj.lb - obj.ub) / 2;
      [obj.miny, id] = min(y);
      obj.minx = X(id, :);
      [obj.maxy, id] = max(y);
      obj.maxx = X(id, :);
    end
    
    function [x, y] = getMax(obj)
      % get point with maximal objective value
      [y, id] = max(obj.y);
      x = obj.X(id, :);
    end
    
    function [x, y] = getMin(obj)
      % get point with minimal objective value
      [y, id] = min(obj.y);
      x = obj.X(id, :);
    end
    
    function d = getDistCtr2Max(obj, distance)
      % get distance from cell center to the point with maximal objective
      % value within the cell
      if nargin < 2
        distance = 'euclidean';
      end
      d = pdist([obj.center; obj.maxx], distance);
    end
    
    function d = getDistCtr2Min(obj, distance)
      % get distance from cell center to the point with minimal objective
      % value within the cell
      if nargin < 2
        distance = 'euclidean';
      end
      d = pdist([obj.center; obj.minx], distance);
    end
    
    function ang = getMaxMinAngle(obj)
      % get angle between the point with maximal objective value, the cell 
      % center, and the point with minimal objective value within the cell
      u = obj.minx - obj.center;
      v = obj.maxx - obj.center;
      ang = acos(dot(u, v) / (norm(u)*norm(v)));
    end
    
    function df = getMaxMinDiff(obj)
      % get difference between the points with minimal and maximal
      % objective values
      df = obj.maxy - obj.miny;
    end
  end
end