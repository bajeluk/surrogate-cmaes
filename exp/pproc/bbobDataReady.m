function data = bbobDataReady(datapath, funcSet)
% Prepares data for further processing.
% Returns cell array 'data' functions x dimensions x settings and 
% appropriate 'settings'.
%
% Input:
%   datapath      - path to data | string
%   funcSet       - structure with fields 'BBfunc' (numbers of BBOB
%                   functions) and 'dims' (numbers of dimensions) 
%                   | structure
%
% Output:
%   data     - aggregated data of size functions x dimensions x settings 
%              | cell array

  if nargin < 2
    if nargin < 1
      data = {};
      help bbobDataReady
      return
    end
    funcSet.BBfunc = 1:24;
    funcSet.dims = 20;
  end

  data = cell(length(funcSet.BBfunc), length(funcSet.dims));
  
  % load results
  assert(isdir(datapath), 'Source path is not a folder')
  
  for f = funcSet.BBfunc
    for d = funcSet.dims
      actualData = importdata(fullfile(datapath, ['data_f', num2str(f)], ['bbobexp_f', num2str(f), '_DIM', num2str(d), '.dat']), ' ', 1);
      
    end
  end

end