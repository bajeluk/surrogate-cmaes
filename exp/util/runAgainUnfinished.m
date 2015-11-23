function runAgainUnfinished(experiment, timebound)
% Function runs again instances of 'experiment' which did not successfully
% ended (according to png file) setting up time boundary 'timebound'

  if nargin < 1
    fprintf('ID of experiment has to be set!\n')
    help runAgainUnfinished
    return
  end
  if nargin < 2
    timebound = '4h';
  end

  expfolder = fullfile('exp', 'experiments', experiment);
  params = load(fullfile(expfolder, 'scmaes_params.mat'));
  nInstances = prod(cellfun(@length, {params.bbParamDef.values, params.sgParamDef.values, params.cmParamDef.values}));
  pngList = dir(fullfile(expfolder, [experiment, '*.png']));
  pngNumbers = sort(cell2mat(cellfun(@(x) str2double(x(strfind(x, 'D_') + 2: strfind(x, '.png') -1)), {pngList(:).name}, 'UniformOutput', false)));
  runAgain = 1:nInstances;
  runAgain(pngNumbers) = [];
  if isempty(runAgain)
    fprintf('All instances have been computed. Nothing to do.\n')
  else
    metacentrum_master_template('exp_restrEC_01', runAgain, timebound);
  end
  
end
