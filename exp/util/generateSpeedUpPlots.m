% script for making speed up graphs on gecco 2015 abstract
% The speed up is considered between the CMA-ES and GP, or RF

exppath = fullfile('exp','experiments');
gppath = fullfile(exppath,'exp_geneEC_06');
rfpath = {[gppath,'_rflite'],[gppath,'_rf5_2'],[gppath,'_rf5_3'],[gppath,'_rf5_4']};

BBfunc = [1,2,3,5,6,8,10,11,12,13,14,20,21];
BBfuncInv = [1,2,3,0,4,5,0,6,0,7,8,9,10,11,0,0,0,0,0,12,13];
dims = [2,5,10];
dimsInv = [0,1,0,0,2,0,0,0,0,3];

rf_evals = cell(13,3);
gp_evals = cell(13,3);

% load and complete results

list = dir(fullfile(gppath,'*.mat'));
gp_list = {list(1:end-1).name};                % get rid of scmaes_params.mat

rf_list = {};
for i = 1:length(rfpath)
  list = dir(fullfile(rfpath{i},'*.mat'));
  if i>1
    rf_list(end+1:end+length(list)-2) = {list(2:end-1).name};
  else
    rf_list(end+1:end+length(list)-1) = {list(1:end-1).name};
  end
end

% gp evaluations load
for i = 1:length(gp_list)
  S = load(gp_list{i},'-mat','y_evals');
  idx = strfind(gp_list{i},'_');
  func = str2double(gp_list{i}(1,idx(end-2)+1:idx(end-1)-1));
  dim = str2double(gp_list{i}(1,idx(end-1)+1:idx(end)-2));
  gp_evals{BBfuncInv(func),dimsInv(dim)} = S.y_evals;
end

% rf evaluations load
for i = 1:length(rf_list)
  idx = strfind(rf_list{i},'_');
  id = str2double(rf_list{i}(1,idx(end)+1:end-4));
  if ~any(strfind(rf_list{i},'rflite')) || (mod(id,6) == 1 && any(strfind(rf_list{i},'rflite')))  % use only data with the same s-cmaes settings
    func = str2double(rf_list{i}(1,idx(end-2)+1:idx(end-1)-1));
    dim = str2double(rf_list{i}(1,idx(end-1)+1:idx(end)-2));
    S = load(rf_list{i},'-mat','y_evals');
    rf_evals{BBfuncInv(func),dimsInv(dim)}(end+1:end+length(S.y_evals),1) = S.y_evals;
  end
end

