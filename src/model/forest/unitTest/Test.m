classdef (Abstract) Test < matlab.unittest.TestCase
  
  properties
    drawEnabled = true
    fig % figure
    iPlot % current plot
    params = struct % parameters of current test
    joinedParams = ''
    name = {'', ''} % name of the current test
  end
  
  methods (TestMethodSetup)
    function setup(testCase)
      rng('default');
      if testCase.drawEnabled
        testCase.fig = figure;
        testCase.iPlot = 0;
      end
    end
  end

  methods (TestMethodTeardown)
    function teardown(testCase)
      if testCase.drawEnabled
        title = sprintf('%s.%s(%s)', ...
          testCase.name{1}, ...
          testCase.name{2}, ...
          testCase.joinedParams);
        set(testCase.fig, 'Name', title, 'NumberTitle', 'off');
        
        path = sprintf('figures/%s', testCase.name{1});
        [~,~,~] = mkdir(path);
        filename = sprintf('%s/%s(%s).png', ...
          path, ...
          testCase.name{2}, ...
          testCase.joinedParams);
        print(testCase.fig, filename, '-dpng');
        close(testCase.fig);
      end
    end
  end
  
  methods (Access = protected)
    function reset(testCase, params, method)
      stack = dbstack;
      name = stack(2).name;
      if nargin >= 3
        name = strcat(name, method);
      end
      testCase.name = strsplit(name, '.');
      testCase.params = params;
      
      testCase.joinedParams = '';
      keys = fieldnames(testCase.params);
      for i = 1:numel(keys)
        value = testCase.params.(keys{i});
        if i == numel(keys)
          format = '%s%s=%s';
        else
          format = '%s%s=%s,';
        end
        testCase.joinedParams = sprintf(format, ...
          testCase.joinedParams, ...
          keys{i}, ...
          value);
      end
    end
  end
  
  methods (Static, Access = protected)
    function drawSplit(X, y, splitIdx)
      for side = [0 1]
        Xs = X(splitIdx == side, :);
        ys = y(splitIdx == side, :);
        [n, d] = size(Xs);
        if n == 0
          continue
        end
        if d == 1
          scatter(Xs(:, 1), ys);
        else
          scatter3(Xs(:, 1), Xs(:, 2), ys);
        end
        hold on;
      end
      hold off;
    end
    
    function [X, minVal, maxVal] = ...
        generateInput(n, d, minVal, maxVal)
      % random points
      if nargin < 1
        n = 1000; % number of points
      end
      if nargin < 2
        d = 1; % dimension
      end
      if nargin < 3
        minVal = -100; % range [minVal, maxVal]
      end
      if nargin < 4
        maxVal = 100; % range [minVal, maxVal]
      end
      X = minVal + (maxVal - minVal) * rand(n, d);
    end
  end
end

