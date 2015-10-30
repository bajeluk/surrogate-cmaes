function runAgainErr(experiment, timebound)
% Function runs again instances of 'experiment' which ended with error
% setting up time boundary 'timebound'

  errorList = dir(fullfile('exp', 'experiments', experiment, '*ERROR.mat'));
  errorNumbers = sort(cell2mat(cellfun(@(x) str2double(x(strfind(x, 'D_') + 2: strfind(x, '_ERROR') -1)), {errorList(:).name}, 'UniformOutput', false)));
  metacentrum_master_template('exp_restrEC_01', errorNumbers, timebound);
  quit(0)
end
