function data = bbobDataReady(datapath, funcSet, maxFE)
% data = bbobDataReady(datapath, funcSet, maxFE) prepares BBOB data for 
% further processing. Returns cell array 'data' functions x dimensions.
%
% Input:
%   datapath - path to data | string
%   funcSet  - structure with fields 'BBfunc' (numbers of BBOB functions) 
%              and 'dims' (numbers of dimensions) | structure
%   maxFE    - maximal number of function evaluations | integer
%
% Output:
%   data     - aggregated data of size functions x dimensions | cell array
%
% See Also:
%   dataReady, divSmooth

  if nargin < 3
    if nargin < 2
      if nargin < 1
        if nargout > 0
          data = {};
        end
        help bbobDataReady
        return
      end
      funcSet.BBfunc = 1:24;
      funcSet.dims = 2;
    end
    maxFE = 250;
  end
  
  if ~iscell(datapath)
    datapath = {datapath};
  end
  nData = length(datapath);

  data = cell(length(funcSet.BBfunc), length(funcSet.dims));
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  
  infoList = [];
  % load results
  [~, workdir] = fileparts(pwd());
  if (~strcmp(workdir, 'surrogate-cmaes'))
    cd([fileparts(mfilename('fullpath')) filesep '..' filesep '..']);
  end
  for dat = 1:nData
    assert(isdir(datapath{dat}), '%s is not a folder', datapath{dat})
    infoList = [infoList; searchFile(datapath{dat}, '*.info')];
  end
  
  % load data from each info file
  for fil = 1:length(infoList)
    actualInfo = importdata(infoList{fil});
    % incomplete data
    if isstruct(actualInfo)
      actualInfo = actualInfo.textdata;    
    end
    infoFileSplit = strsplit(infoList{fil}, filesep);
    if isempty(infoFileSplit{1})
      infoFileSplit{1} = filesep;
    end

    dataRowID = find(cellfun(@(x) ~isempty(strfind(x, '.dat')), actualInfo(:, 1)));
    % load .dat files
    for r = 1:length(dataRowID)
      % complete row
      if all(cellfun(@isempty, actualInfo(dataRowID(r), 2:end)))
        rowSplit = strsplit(actualInfo{dataRowID(r)}, ', ');
      % incomplete row
      else
        rowSplit = actualInfo(dataRowID(r), :);
      end
      % name of data file
      datName = strrep(strrep(rowSplit{1}, '\', filesep), '/', filesep);
      datFile = fullfile(infoFileSplit{1:end-1}, datName);
      tdatFile = [datFile(1:end-3), 'tdat'];
      if exist(tdatFile, 'file')
      % uncomment for .tdat file data loading - has a bug somewhere
%         datFile = tdatFile;
      end
      % extract function and dimension number
      datSplit = strsplit(rowSplit{1}, '_');
      f = str2double(datSplit{end-1}(2:end));
      d = str2double(datSplit{end}(4:end-4));
      % numbers of instances
      instanceNum = cellfun(@(x) str2double(x(1 : strfind(x, ':')-1)), rowSplit(2:end));
      if exist(datFile, 'file') && any(funcSet.BBfunc == f) && any(funcSet.dims == d)
        fi = funcSet.BBfuncInv(f);
        di = funcSet.dimsInv(d);
        actualData = importdata(datFile, ' ', 1);
        % find individual instances
        % TODO: find instances properly - now supposing each instance starts
        % with 1 evaluation
        instanceStartId = inverseIndex(actualData.data(:, 1) == 1);
        for i = 1 : length(instanceStartId)-1
          data{fi, di}{instanceNum(i)} = actualData.data(instanceStartId(i):instanceStartId(i+1)-1, [3,1]);
        end
        data{fi, di}{instanceNum(end)} = actualData.data(instanceStartId(end):end, [3,1]);
      end
    end
  end
    
  % smooth data
  data = divSmooth(data, funcSet, maxFE);
  
end