function data = readSMACResults(datapath, funcSet)
% Returns SMAC results in cell array 'data' functions x dimensions x 
% settings and appropriate 'settings'.
%
% Input:
%   datapath - path to data | string
%   funcSet  - structure with fields 
%                'BBfunc' (numbers of BBOB functions) 
%                'dims' (numbers of dimensions)
%                'instances' (numbers of function instances)
%
% Output:
%   data     - aggregated data of size functions x dimensions x settings 
%              | cell array

  if nargin < 2
    if nargin < 1
      data = {};
      help readSMACResults
      return
    end
    funcSet.BBfunc = 1:24;
    funcSet.dims = 2;
  end

  data = cell(length(funcSet.BBfunc), length(funcSet.dims));
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  funcSet.instances = defopts(funcSet, 'instances', [1:5, 31:40]);
  
  % load results
  assert(isdir(datapath), 'Source path is not a folder')
  
  for f = funcSet.BBfunc
    for d = funcSet.dims
      for i = 1:length(funcSet.instances)
        try
          actualData = importdata(fullfile(datapath, ...
            ['dim', num2str(d), '-factor100iinstance', num2str(funcSet.instances(i))], ...
            ['data_f', num2str(f)], ...
            ['bbobexp_f', num2str(f), '_DIM', num2str(d), '.dat']), ' ', 1);
          fi = funcSet.BBfuncInv(f);
          di = funcSet.dimsInv(d);
          data{fi, di}{i} = actualData.data(:, [3,1]);
        catch
          warning('Unable to open function %d dimension %d', f, d)
        end
      end
    end
  end
  
  data = divSmooth(data, funcSet);

end