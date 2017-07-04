%FILLTMPRESULTS -- re-save old version of *_tmp_* results in *_results_* format
%
% fillTmpResults(exp_id, exppath_short)
%
% Traverses through directory [exppath_short '/' exp_id'], loads each file
%   ${EXP_ID}_tmp_[0-9]+\.mat
% and saves into 
%   ${EXP_ID}_tmp_${FUN}_${DIM}D_[0-9]+\.mat
% also the bbParams, surrogateParams, cmaesParams structures
function fillTmpResults(exp_id, exppath_short)
  exppath = [exppath_short filesep exp_id]; 
  cd(exppath);
  load([exppath filesep 'scmaes_params.mat']);

  thisDir = dir('*_tmp_*.mat');

  for i = 1:length(thisDir)
    if (isempty(regexp(thisDir(i).name, [exp_id '_tmp_[0-9]+\.mat'])))
      continue;
    else
      id = str2num(regexprep(thisDir(i).name, [exp_id '_tmp_([0-9]+)\.mat'], '$1'));
    end

    disp(['loading ' thisDir(i).name '... (' num2str(id) ')']);

    [bbParams, surrogateParams, cmaesParams, nNonBbobValues] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

    load(thisDir(i).name);

    expFileID = [num2str(exp_settings.bbob_function) '_' num2str(exp_settings.dim) 'D_' num2str(id)];
    newTmpFile = [exppath filesep exp_settings.exp_id '_tmp_' expFileID '.mat'];
    disp(['saving into ' newTmpFile ': exp_id, exp_settings, exp_results, y_evals, surrogateParams, cmaesParams, bbParams, cmaes_out']);
    save(newTmpFile, 'exp_id', 'exp_settings', 'exp_results', 'y_evals', 'surrogateParams', 'cmaesParams', 'bbParams', 'cmaes_out');
  end
end
