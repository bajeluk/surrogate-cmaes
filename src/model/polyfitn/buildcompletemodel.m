function [modelterms,varlist] = buildcompletemodel(order,p)
% 
% arguments: (input)
%  order - scalar integer, defines the total (maximum) order 
%
%  p     - scalar integer - defines the dimension of the
%          independent variable space
%
% arguments: (output)
%  modelterms - exponent array for the model
%
%  varlist - cell array of character variable names

% build the exponent array recursively
if p == 0
  % terminal case
  modelterms = [];
elseif (order == 0)
  % terminal case
  modelterms = zeros(1,p);
elseif (p==1)
  % terminal case
  modelterms = (order:-1:0)';
else
  % general recursive case
  modelterms = zeros(0,p);
  for k = order:-1:0
    t = buildcompletemodel(order-k,p-1);
    nt = size(t,1);
    modelterms = [modelterms;[repmat(k,nt,1),t]];
  end
end

% create a list of variable names for the variables on the fly
varlist = cell(1,p);
for i = 1:p
  varlist{i} = ['X',num2str(i)];
end
