function data = bbobDataReady(datapath, funcSet)
% Prepares data for further processing.
% Returns cell array 'data' functions x dimensions x settings and 
% appropriate 'settings'.
%
% Input:
%   datapath - path to data | string
%   funcSet  - structure with fields 'BBfunc' (numbers of BBOB functions) 
%              and 'dims' (numbers of dimensions) | structure
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
    funcSet.dims = 2;
  end

  data = cell(length(funcSet.BBfunc), length(funcSet.dims));
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  
  % load results
  assert(isdir(datapath), 'Source path is not a folder')
  
  for f = funcSet.BBfunc
    for d = funcSet.dims
      try
        actualData = importdata(fullfile(datapath, ['data_f', num2str(f)], ['bbobexp_f', num2str(f), '_DIM', num2str(d), '.dat']), ' ', 1);
        % find individual instances
        % TODO: find instances properly - now supposing each instance starts
        % with 1 evaluation
        instanceStartId = inverseIndex(actualData.data(:, 1) == 1);
        fi = funcSet.BBfuncInv(f);
        di = funcSet.dimsInv(d);
        for i = 1 : length(instanceStartId)-1
          data{fi, di}{i} = actualData.data(instanceStartId(i):instanceStartId(i+1)-1, [3,1]);
        end
        data{fi, di}{length(instanceStartId)} = actualData.data(instanceStartId(end):end, [3,1]);
      catch
        warning('Unable to open function %d dimension %d', f, d)
      end
    end
  end
  
  data = divSmooth(data, funcSet);
  
end