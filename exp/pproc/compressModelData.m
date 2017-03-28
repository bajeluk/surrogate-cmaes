function compressedFolder = compressModelData(modelFolder, modelPreposition)
% compressModelData - removes models from results of model testing,
%                     i.e. modelPreposition*.mat files saved in
%                     modelFolder, and saves them into modelFolder_compressed
%
% Input:
%   modelFolder - folder with files from testModels | string
%   modelPreposition - only files whose names start with this string will be
%                      compressed | string
% Output:
%   modelFolder  - folder containing results | string
%
% based on compressModelData() by Jan Juranko, 2017
%          https://github.com/juranja3/surrogate-cmaes/blob/modeltest/exp/compressModelData.m

  if (~exist('modelPreposition', 'var'))
    modelPreposition = [];
  end

  compressedFolder = [modelFolder, '_compressed'];
  fprintf('Compressing data in ''%s'', output folder is ''%s''\n', modelFolder, compressedFolder);

  if (exist(compressedFolder, 'dir'))
    warning('Output directory already exist! You have 5 sec to stop rewriting...');
    pause(5);
  else
    mkdir(compressedFolder);
  end

  % data are saved in separate folders
  dirs = dir(fullfile(modelFolder, '*model_*'));
  dirIndex = find([dirs.isdir]);

  for i = 1:length(dirIndex)
    dirName = dirs(dirIndex(i)).name;
    fprintf('Processing ''%s''\n', dirName);
    mo_struct = [];

    modelFiles = dir(fullfile(modelFolder, dirName, [modelPreposition, '*.mat']));
    % for each file that matches create its compressed copy in compressedFolder
    for j = 1:length(modelFiles)
      fprintf('... file ''%s''\n', modelFiles(j).name);
      outFile = fullfile(compressedFolder,dirName, modelFiles(j).name);
      if (exist(outFile, 'file'))
        fprintf('... ... Processing skipped, the compressed file already exists.\n');
        continue;
      end

      % load original data
      data = load(fullfile(modelFolder,dirName, modelFiles(j).name));

      % identify the right modelOption entry during the first file processing
      if (isempty(mo_struct))
        mo_struct = getThisModelOption(dirName, data.modelOptions);
      end

      % compress (delete dataset in current implementation)
      % data = compressData(data);

      % omit the models totally
      data = rmfield(data, 'models');

      % save only *current* modelOptions struc (not all of them)
      data.modelOptions = mo_struct;

      % create folder
      [~,~] = mkdir(fullfile(compressedFolder,dirName));
      % save compressed data
      save(outFile, '-struct', 'data');
    end
  end
end

function data = compressData(data)
  data.modelOptions = [];
  for i=size(data.models,1)
    for j=size(data.models,2)
      data.models{i,j}.dataset=[];
    end
  end
end

function [mo_struct, mo_idx] = getThisModelOption(dirName, modelOptions)
  m = regexp(dirName, '_([0-9]+)_', 'tokens');
  if (~isempty(m{1}{1}))
    dirHash = str2double(m{1}{1});

    for i = 1:length(modelOptions)
      if (str2double(modelHash(modelOptions{i})) == dirHash)
        mo_idx = i;
        mo_struct = modelOptions{i};
        return;
      end
    end
  end
  fprintf(2, 'The right index of the directory hash not found in modelOptions.\n');
  mo_idx = -1;
  mo_struct = [];
end
