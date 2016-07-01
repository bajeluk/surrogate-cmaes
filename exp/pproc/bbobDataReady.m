function data = bbobDataReady(datapath, funcSet)
% Prepares data for further processing.
% Returns cell array 'data' functions x dimensions x settings and 
% appropriate 'settings'.
% (Warning: Suppose that it does not matter how instances are ordered.)
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
  
  if ~iscell(datapath)
    datapath = {datapath};
  end
  nData = length(datapath);

  data = cell(length(funcSet.BBfunc), length(funcSet.dims));
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  
  % load results
  for dat = 1:nData
    assert(isdir(datapath{dat}), '%s is not a folder', datapath{dat})
  end
  
  for f = funcSet.BBfunc
    for d = funcSet.dims
      usefulFiles = false(1, nData);
      filename = cell(1, nData);
      for dat = 1:nData
        filename{dat} = fullfile(datapath{dat}, ['data_f', num2str(f)], ['bbobexp_f', num2str(f), '_DIM', num2str(d), '.dat']);
        if exist(filename{dat}, 'file')
          usefulFiles(dat) = true;
        end
      end
      if any(usefulFiles)
        fi = funcSet.BBfuncInv(f);
        di = funcSet.dimsInv(d);
        % supposing that it does not matter how instances are ordered
        for dat = inverseIndex(usefulFiles)
          actualData = importdata(filename{dat}, ' ', 1);
          % find individual instances
          % TODO: find instances properly - now supposing each instance starts
          % with 1 evaluation
          instanceStartId = inverseIndex(actualData.data(:, 1) == 1);
          actualLength = length(data{fi, di});
          for i = 1 : length(instanceStartId)-1
            data{fi, di}{actualLength + i} = actualData.data(instanceStartId(i):instanceStartId(i+1)-1, [3,1]);
          end
          data{fi, di}{actualLength + length(instanceStartId)} = actualData.data(instanceStartId(end):end, [3,1]);
        end
      else
        warning('Unable to open function %d dimension %d.', f, d)
      end
    end
  end
  
  data = divSmooth(data, funcSet);
  
end