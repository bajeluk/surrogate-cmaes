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
    function obj = CMCell(X, y, dim, lb, ub)
      % constructor
        if nargin < 1
          X = [];
          y = [];
          dim = 0;
          lb = [];
          ub = [];
        end
      % empty set
      if isempty(X)
        obj.X = [];
        obj.y = [];
        obj.dim = dim;
        obj.lb = lb;
        obj.ub = ub;
        obj.miny = [];
        obj.minx = [];
        obj.maxy = [];
        obj.maxx = [];
      else
        obj.X = X;
        obj.y = y;
        obj.dim = dim;
        obj.lb = lb;
        obj.ub = ub;
        [obj.miny, id] = min(y);
        obj.minx = X(id, :);
        [obj.maxy, id] = max(y);
        obj.maxx = X(id, :);
      end
      obj.center = obj.lb + (obj.lb - obj.ub) / 2;
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
    
    function y_mean = getMean(obj)
      % get mean objective value
      y_mean = mean(obj.y);
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
      if isempty(obj.X)
        ang = [];
      else
        u = obj.minx - obj.center;
        v = obj.maxx - obj.center;
        ang = abs(acos(dot(u, v) / (norm(u)*norm(v))));
      end
    end
    
    function df = getMaxMinDiff(obj)
      % get difference between the points with minimal and maximal
      % objective values
      if isempty(obj.X)
        df = [];
      else
        df = obj.maxy - obj.miny;
      end
    end
    
    function [x, y] = getNearCtrPoint(obj, cl_distance, dist_param)
      % get point nearest to the cell center
      if isempty(obj.X)
        x = [];
        y = [];
      else
        % default distance
        if nargin < 2
          cl_distance = 'euclidean';
        end
        % minkowski and mahalanobis settings
        if any(strcmp(cl_distance, {'minkowski', 'mahalanobis'})) && nargin == 3
          [~, id] = pdist2(obj.X, obj.center, cl_distance, dist_param, 'Smallest', 1);
        % other distances
        else
          [~, id] = pdist2(obj.X, obj.center, cl_distance, 'Smallest', 1);
        end
        x = obj.X(id, :);
        y = obj.y(id);
      end
    end
    
    function gradHomo = getGradHomogeneity(obj, cl_distance, dist_param)
      % get gradient homogeneity of the cell
      nPoints = numel(obj.y);
      
      if nPoints < 3
        % 2 or less points are not helpful
        gradHomo = [];
      else
        % default distance
        if nargin < 2
          cl_distance = 'euclidean';
        end
        % minkowski and mahalanobis settings
        if any(strcmp(cl_distance, {'minkowski', 'mahalanobis'})) && ...
            nargin == 3 && ~isempty(dist_param)
          [pDistances, id] = pdist2(obj.X, obj.X, cl_distance, dist_param, 'Smallest', 2);
        % other distances
        else
          [pDistances, id] = pdist2(obj.X, obj.X, cl_distance, 'Smallest', 2);
        end
        % the second row contains required ids, the first row are point ids
        % themselves
        pDistances = pDistances(2, :)';
        % calculate "gradient"
        y_direction = (obj.y(id(1,:)) < obj.y(id(2, :)));
        % calculate direction
        X_direction = obj.X(id(2,:), :) - obj.X(id(1,:), :);
        % compute normalized vectors
        nv = repmat(((1./pDistances) .* (2*y_direction - 1)), 1, obj.dim) .* X_direction;
        % calculate lenght of sum of all vectors
        gradHomo = norm(sum(nv)) / nPoints;
      end
    end
    
    function lm = fitPolyModel(obj, modelspec)
      % fit polynomial model on cell points
      if nargin < 2
        modelspec = 'linear';
      end
      if numel(obj.y) >= getMinFitlmPoints(modelspec, obj.dim) && ...
          ~all(isnan(obj.y))
        lm = fitlm(obj.X, obj.y, modelspec);
      else
        % empty linear model
        lm = LinearModel();
      end
    end
    
    function state = isempty(obj)
      % cell is empty when it does not contain any data
      state = isempty(obj.X);
    end
    
  end
end

function N = getMinFitlmPoints(modelspec, dim)
% returns minimal number of points necessary to train the model through
% fitlm function
  switch modelspec
    case 'constant'
      % could be only one, but fitlm fails having only one point
      N = 2;
    case 'linear'
      N = 1 + dim;
    case 'interactions'
      N = 1 + dim + dim.*(dim-1)/2;
    case 'purequadratic'
      N = 1 + 2*dim;
    case 'quadratic'
      N = 1 + 2*dim + dim.*(dim-1)/2;
    otherwise
      % poly type
      if strcmp(modelspec(1:4), 'poly')
        % gain degrees
        deg = cell(numel(modelspec) - 4, 1);
        for i = 1 : numel(modelspec) - 4
          deg{i} = 0 : str2double(modelspec(i - 4));
        end
        % find all combinations
        combinations = combvec(deg{:});
        % count only combinations with degree lower than given maximal 
        % degree
        N = sum(sum(combinations) <= max(max(combinations)));
      else
        error('Cannot interpret ''%s'' as model formula.', modelspec)
      end
  end
end