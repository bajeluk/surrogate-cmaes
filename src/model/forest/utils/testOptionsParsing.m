function t = testOptionsParsing()
    args = {'option1', 'option2'};
    fcns = {
        @foo_isfield ;
        @foo_inputparser
    };

    % parameters sweep
    [f,k] = ndgrid(1:numel(fcns), 0:numel(args));
    f = f(:); k = k(:);

    % test combinations of functions and number of input args
    t = cell(numel(f), 3);
    for i=1:size(t,1)
        t{i,1} = func2str(fcns{f(i)});
        t{i,2} = k(i);
        if k(i) == 0
          t{i,3} = timeit(@() feval(fcns{f(i)}), 2);
        elseif k(i) == 1
          t{i,3} = timeit(@() feval(fcns{f(i)}, struct('option1', 1)), 2);
        elseif k(i) == 2
          t{i,3} = timeit(@() feval(fcns{f(i)}, struct('option1', 1, 'option2', 2)), 2);
        end
    end

    % format results in table
    t = cell2table(t, 'VariableNames',{'func','nargs','time'});

    figure(2);
    tt = unstack(t, 'time', 'func');
    names = tt.Properties.VariableNames(2:end);
    bar(tt{:,2:end}.')
    set(gca, 'XTick',1:numel(names), 'XTickLabel',names, 'YGrid','on')
    legend(num2str(tt{:,1}, 'nargin=%d'))
    ylabel('Time [sec]'), xlabel('Functions')
end

function [aa,bb] = foo_isfield(options)
  if nargin < 1, options = struct(); end
  if ~isfield(options, 'option1'); options.option1 = 1; end;
  if ~isfield(options, 'option2'); options.option2 = 1; end;
  aa = options.option1;
  bb = options.option2;
end

function [aa,bb] = foo_inputparser(varargin)
  p = inputParser;
  p.addParamValue('option1', 1);
  p.addParamValue('option2', 1);
  p.parse(varargin{:});
  aa = p.Results.option1;
  bb = p.Results.option2;
end