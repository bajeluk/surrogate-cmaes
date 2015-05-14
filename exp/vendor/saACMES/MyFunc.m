function Fit = MyFunc(x)
global settings;

Fit = feval('fgeneric',x) - settings.ftarget;    % does not change the results, but better for logs
%disp(x);
