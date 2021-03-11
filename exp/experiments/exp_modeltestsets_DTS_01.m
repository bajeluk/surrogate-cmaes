% script creating testing sets for surrogate models from the runs of the
% DTS-CMA-ES in exp_doubleEC_28_log_nonadapt

% set seed for replicable results

% basic settings
exp_id = 'exp_doubleEC_28_log_nonadapt';
fun = 1:24;
dim = [2, 3, 5, 10, 20];
inst = 11:15;
% name settings
datasetName = 'DTS_meta_005_new';
outputDirname = exp_id;
% process settings
isForData = true;
loadModels = false;
nSnapshotsPerRun = 25;
rewriteResults = false;
sampleMethod = 'uniform_wor'; % uniform sampling without replacement

% create testing sets
ds = modelTestSets(exp_id, fun, dim, inst, ...
                   'datasetName', datasetName, ...
                   'isForData', isForData, ...
                   'loadModels', loadModels, ...
                   'nSnapshotsPerRun', nSnapshotsPerRun, ...
                   'outputDirname', outputDirname, ...
                   'rewriteResults', rewriteResults, ...
                   'sampleMethod', sampleMethod ...
                  );