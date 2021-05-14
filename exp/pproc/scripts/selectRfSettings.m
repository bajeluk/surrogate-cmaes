% Script for selection of the most convenient RF settings for SMS problem
% investigation.
% RF models differ in decision tree splitting method settings (CART, OC1, 
% PAIR, SECRET, SUPPORT), number of trees in RF (64, 128, 256, 512, 1024),
% number of tree points ({0.25, 0.5, 0.75, 1}*number of available points), 
% and number of tree dimensions ({0.25, 0.5, 0.75, 1}*number of available 
% dimensions).
% 100 RF settings selected out of possible 400 using lhs design was tested
% on DTS_meta_005_validation.mat utilizing TSS nearest and TSS full.
% Using results from exp_DTSmodels_meta_03_rf_validation_lhs and
% exp_DTSmodels_meta_03_rf_validation_full_lhs.

%%

% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'selectRfSettings_data.mat');
if isfile(tmpFName)
  load(tmpFName);
else
  
%%

% properties
dimensions = [2, 3, 5, 10, 20];
functions = 1:24;
instances = 11:15;
generations = 1:25; % generation snapshots
splitFunc = { @AxisSplit, @GaussianSplit, @HillClimbingObliqueSplit, @PairObliqueSplit, @ResidualObliqueSplit };
nTrees = 2.^[6, 7, 8, 9, 10];
nFeaturesToSample = { 'ceil(0.25*dim)', 'ceil(0.5*dim)', 'ceil(0.75*dim)', 'dim' };
inBagFraction = [0.25, 0.5, 0.75, 1.0];

% experiment paths
expname = 'exp_DTSmodels_meta_03_rf_validation_lhs';
expname_full = 'exp_DTSmodels_meta_03_rf_validation_full_lhs';
expfolder = fullfile('exp', 'experiments', expname);
expfolder_full = fullfile('exp', 'experiments', expname_full);
resfolder = fullfile('exp', 'experiments', expname);
modelDirList = dir(fullfile(expfolder, 'forestmodel_*'));
modelDirList_full = dir(fullfile(expfolder_full, 'forestmodel_*'));

% init results
resSize = [numel(dimensions), numel(functions), numel(instances), ...
           numel(generations), numel(splitFunc), numel(nTrees), ...
           numel(nFeaturesToSample), numel(inBagFraction)];
mseRes = NaN(resSize);
kendallRes = NaN(resSize);
rdeRes = NaN(resSize);
mseRes_full = NaN(resSize);
kendallRes_full = NaN(resSize);
rdeRes_full = NaN(resSize);

coorCounter = zeros(1, prod(resSize));
coorCounter_full = zeros(1, prod(resSize));
fileNum = 0;

% TSS nearest model settings loop
for m = 1:numel(modelDirList)
  fileDirList = dir(fullfile(expfolder, modelDirList(m).name, 'forestmodel_*.mat'));
  % file loop
  for fil = 1:numel(fileDirList)
    tic;
    % load file
    S = load(fullfile(fileDirList(fil).folder, fileDirList(fil).name));
    fileNum = fileNum + 1;
    % inst loop
    for inst = 1:size(S.stats.mse, 1)
      % generation snapshot loop
      for g = 1:size(S.stats.mse, 3)
        % add function, dimension, instance, generation snapshot, splitting
        % function, number of trees, number of features to sample, and
        % number of bagged points
        resCoor = [find(S.dim == dimensions), ...
                   find(S.fun == functions), ...
                   find(S.instances(inst) == instances), ...
                   find(g == generations), ...
                   find(cellfun(@(x) isequal(S.modelOptions1.tree_splitFunc, x), splitFunc)), ...
                   find(S.modelOptions1.rf_nTrees == nTrees), ...
                   find(strcmp(S.modelOptions1.rf_nFeaturesToSample, nFeaturesToSample)), ...
                   find(S.modelOptions1.rf_inBagFraction == inBagFraction)];
%         if ~isempty(coorList) && ismember(resCoor, coorList, 'rows')
%           return
%         end
%         coorList = [coorList; resCoor];
        resInd = coor2ind(resCoor, resSize);
        coorCounter(resInd) = coorCounter(resInd) + 1;
        % add error
        mseRes(resInd) = S.stats.mse(inst, 1, g);
        kendallRes(resInd) = S.stats.kendall(inst, 1, g);
        rdeRes(resInd) = S.stats.rde(inst, 1, g);
      end
    end
    t = toc;
    % print message
    fprintf('%6d/12000: %0.2fs %s\n', fileNum, t, fullfile(fileDirList(fil).folder, fileDirList(fil).name))
  end
end

%% TSS full model settings loop
fileNum = 0;
for m = 1:numel(modelDirList_full)
  fileDirList = dir(fullfile(expfolder_full, modelDirList_full(m).name, 'forestmodel_*.mat'));
  % file loop
  for fil = 1:numel(fileDirList)
    tic;
    % load file
    S = load(fullfile(fileDirList(fil).folder, fileDirList(fil).name));
    fileNum = fileNum + 1;
    % inst loop
    for inst = 1:size(S.stats.mse, 1)
      % generation snapshot loop
      for g = 1:size(S.stats.mse, 3)
        % add function, dimension, instance, generation snapshot, splitting
        % function, number of trees, number of features to sample, and
        % number of bagged points
        resCoor = [find(S.dim == dimensions), ...
                   find(S.fun == functions), ...
                   find(S.instances(inst) == instances), ...
                   find(g == generations), ...
                   find(cellfun(@(x) isequal(S.modelOptions1.tree_splitFunc, x), splitFunc)), ...
                   find(S.modelOptions1.rf_nTrees == nTrees), ...
                   find(strcmp(S.modelOptions1.rf_nFeaturesToSample, nFeaturesToSample)), ...
                   find(S.modelOptions1.rf_inBagFraction == inBagFraction)];
%         if ~isempty(coorList) && ismember(resCoor, coorList, 'rows')
%           return
%         end
%         coorList = [coorList; resCoor];
        resInd = coor2ind(resCoor, resSize);
        coorCounter_full(resInd) = coorCounter_full(resInd) + 1;
        % add error
        mseRes_full(resInd) = S.stats.mse(inst, 1, g);
        kendallRes_full(resInd) = S.stats.kendall(inst, 1, g);
        rdeRes_full(resInd) = S.stats.rde(inst, 1, g);
      end
    end
    t = toc;
    % print message
    fprintf('%6d/12000: %0.2fs %s\n', fileNum, t, fullfile(fileDirList(fil).folder, fileDirList(fil).name))
  end
end

% clear unnecessary variables
clear fil fileDirList fileNum g inst m OCT resTable S t

% save loaded data
save(tmpFName)

end

%% Get best settings for each split method in TSS nearest

coorInd = find(coorCounter);
coorList = ind2coor(coorInd, resSize);
[uniqSet, ~, coorId] = unique(coorList(:, 5:8), 'rows');

% settings loop
for s = 1:max(coorId)
  actCoor = coor2ind(coorList(s==coorId, :), resSize);
  % calculate rde
  rdeMeans(s) = nanmean(rdeRes(actCoor));
  rdeNans(s) = sum(isnan(rdeRes(actCoor))) / numel(rdeRes(actCoor));
  % calculate mse
  mseMeans(s) = nanmean(mseRes(actCoor));
  mseNans(s) = sum(isnan(mseRes(actCoor))) / numel(mseRes(actCoor));
end

clear actCoor s
save(tmpFName)

%% Get best settings for each split method in TSS full

coorInd_full = find(coorCounter_full);
coorList_full = ind2coor(coorInd_full, resSize);
[uniqSet_full, ~, coorId_full] = unique(coorList_full(:, 5:8), 'rows');

% settings loop
for s = 1:max(coorId_full)
  actCoor = coor2ind(coorList_full(s==coorId_full, :), resSize);
  % calculate rde
  rdeMeans_full(s) = nanmean(rdeRes_full(actCoor));
  rdeNans_full(s) = sum(isnan(rdeRes_full(actCoor))) / numel(rdeRes_full(actCoor));
  % calculate mse
  mseMeans_full(s) = nanmean(mseRes_full(actCoor));
  mseNans_full(s) = sum(isnan(mseRes_full(actCoor))) / numel(mseRes_full(actCoor));
end

clear actCoor s
save(tmpFName)

%% Plot resulting errors and numbers of NaNs
dispChar = {'x'};
methodColor = {'m', 'g', 'b', 'y', 'c'};
methodName = {'Axis', 'Gauss', 'HillClimb', 'Pair', 'Residual'};
nFeatToSam = {'D/4', 'D/2', '3D/4', 'D'};
nanTreshold = 0.25;

% TSS nearest RDE scatter
figure()
hold on
% split method loop
for sm = 1:numel(splitFunc)
  calcRange = (sm-1)*20+1:sm*20;
  % plot all settings
  h(sm) = scatter(rdeMeans(calcRange), rdeNans(calcRange), 30, methodColor{sm}, 'filled');
  % mark minimal setting
  penalizedError = rdeMeans(calcRange) + (rdeNans(calcRange) > nanTreshold);
  [~, minSetIdRde(sm)] = min(penalizedError);
  minSetIdRde(sm) = minSetIdRde(sm) + (sm-1)*20;
  minPointCoor = [rdeMeans(minSetIdRde(sm)), rdeNans(minSetIdRde(sm))];
  scatter(minPointCoor(1), minPointCoor(2), 80, 'k', 'filled')
  scatter(minPointCoor(1), minPointCoor(2), 30, methodColor{sm}, 'filled')
  % put best values to description
  descText = ['[', ...
              num2str(nTrees(uniqSet(minSetIdRde(sm), 2))), ',', ...
              nFeatToSam{uniqSet(minSetIdRde(sm), 3)}, ',', ...
              num2str(inBagFraction(uniqSet(minSetIdRde(sm), 4))), ']'];
%   text(minPointCoor(1), minPointCoor(2) - 0.002, descText)
  methodLeg{sm} = [methodName{sm}, ' ', descText];
end
% add treshold line
line(get(gca, 'XLim')+[eps, -eps], [nanTreshold, nanTreshold], 'Color', 'red')
title('RF settings TSS nearest')
xlabel('RDE')
ylabel('% of NaNs')
legend(h, methodLeg, 'Location', 'northwest')
hold off

% TSS nearest MSE scatter
figure()
hold on
% split method loop
for sm = 1:numel(splitFunc)
  calcRange = (sm-1)*20+1:sm*20;
  % plot all settings
  h(sm) = scatter(mseMeans(calcRange), mseNans(calcRange), 30, methodColor{sm}, 'filled');
  % mark minimal setting
  penalizedError = mseMeans(calcRange) + max(mseMeans)*(mseNans(calcRange) > nanTreshold);
  [~, minSetIdMse(sm)] = min(penalizedError);
  minSetIdMse(sm) = minSetIdMse(sm) + (sm-1)*20;
  minPointCoor = [mseMeans(minSetIdMse(sm)), mseNans(minSetIdMse(sm))];
  scatter(minPointCoor(1), minPointCoor(2), 80, 'k', 'filled')
  scatter(minPointCoor(1), minPointCoor(2), 30, methodColor{sm}, 'filled')
  descText = ['[', ...
              num2str(nTrees(uniqSet(minSetIdMse(sm), 2))), ',', ...
              nFeatToSam{uniqSet(minSetIdMse(sm), 3)}, ',', ...
              num2str(inBagFraction(uniqSet(minSetIdMse(sm), 4))), ']'];
%   text(minPointCoor(1), minPointCoor(2) - 0.002, descText)
  methodLeg{sm} = [methodName{sm}, ' ', descText];
end
% add treshold line
line(get(gca, 'XLim')+[eps, -eps], [nanTreshold, nanTreshold], 'Color', 'red')
title('RF settings TSS nearest')
xlabel('MSE')
ylabel('% of NaNs')
legend(h, methodLeg, 'Location', 'northwest')
hold off

% TSS full RDE scatter
figure()
hold on
% split method loop
for sm = 1:numel(splitFunc)
  calcRange = (sm-1)*20+1:sm*20;
  % plot all settings
  h(sm) = scatter(rdeMeans_full(calcRange), rdeNans_full(calcRange), 30, methodColor{sm}, 'filled');
  % mark minimal setting
  penalizedError = rdeMeans_full(calcRange) + (rdeNans_full(calcRange) > nanTreshold);
  [~, minSetIdRde_full(sm)] = min(penalizedError);
  minSetIdRde_full(sm) = minSetIdRde_full(sm) + (sm-1)*20;
  minPointCoor = [rdeMeans_full(minSetIdRde_full(sm)), rdeNans_full(minSetIdRde_full(sm))];
  scatter(minPointCoor(1), minPointCoor(2), 80, 'k', 'filled')
  scatter(minPointCoor(1), minPointCoor(2), 30, methodColor{sm}, 'filled')
  descText = ['[', ...
              num2str(nTrees(uniqSet(minSetIdRde_full(sm), 2))), ',', ...
              nFeatToSam{uniqSet(minSetIdRde_full(sm), 3)}, ',', ...
              num2str(inBagFraction(uniqSet(minSetIdRde_full(sm), 4))), ']'];
%   text(minPointCoor(1), minPointCoor(2) - 0.002, descText)
  methodLeg{sm} = [methodName{sm}, ' ', descText];
end
% add treshold line
line(get(gca, 'XLim')+[eps, -eps], [nanTreshold, nanTreshold], 'Color', 'red')
title('RF settings TSS full')
xlabel('RDE')
ylabel('% of NaNs')
legend(h, methodLeg, 'Location', 'northeast')
hold off

% TSS full MSE scatter
figure()
hold on
% split method loop
for sm = 1:numel(splitFunc)
  calcRange = (sm-1)*20+1:sm*20;
  % plot all settings
  h(sm) = scatter(mseMeans_full(calcRange), mseNans_full(calcRange), 30, methodColor{sm}, 'filled');
  % mark minimal setting
  penalizedError = mseMeans_full(calcRange) + max(mseMeans)*(mseNans_full(calcRange) > nanTreshold);
  [~, minSetIdMse_full(sm)] = min(penalizedError);
  minSetIdMse_full(sm) = minSetIdMse_full(sm) + (sm-1)*20;
  minPointCoor = [mseMeans_full(minSetIdMse_full(sm)), mseNans_full(minSetIdMse_full(sm))];
  scatter(minPointCoor(1), minPointCoor(2), 80, 'k', 'filled')
  scatter(minPointCoor(1), minPointCoor(2), 30, methodColor{sm}, 'filled')
  descText = ['[', ...
              num2str(nTrees(uniqSet(minSetIdMse_full(sm), 2))), ',', ...
              nFeatToSam{uniqSet(minSetIdMse_full(sm), 3)}, ',', ...
              num2str(inBagFraction(uniqSet(minSetIdMse_full(sm), 4))), ']'];
%   text(minPointCoor(1), minPointCoor(2) - 0.002, descText)
  methodLeg{sm} = [methodName{sm}, ' ', descText];
end
% add treshold line
line(get(gca, 'XLim')+[eps, -eps], [nanTreshold, nanTreshold], 'Color', 'red')
title('RF settings TSS full')
xlabel('MSE')
ylabel('% of NaNs')
legend(h, methodLeg, 'Location', 'northwest')
hold off

%% 

clear calcRange descText methodLeg minPointCoor penalizedError sm
% save results
save(fullfile(resfolder, 'selectedRfSettings.mat'))