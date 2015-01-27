function [resRmse, resCorr] = modelTrainTest(modelStr, modelOpts, filelist)
  % run training and testing of model-regression
  %
  % modelStr  -- string of the model to be tested        'gp' | 'rf'
  % modelOpts -- struct with model-parameters
  % filelist  -- string with the filename with list of files
  %              to be loaded with train and testing sets
  %
  % returns:
  %   resRmse -- 3D array with RMSE on the testing sets
  %   resCorr -- 3D array with Kendall's correlations on the testing sets
  
  % load list of training mat-files
  fid = fopen(filelist, 'r');
  files = cell(0);
  fnums = [];
  dims = [];
  evals = [];
  i = 1;
  tline = fgetl(fid);
  while (ischar(tline))
    if (isempty(regexp(tline, '^#', 'ONCE')))
      files{i} = tline;
      % collect function numbers and dimensions
      [regstrings] = regexpi(tline, '_f(\d+)_(\d+)D_\d+_(\d+)\.mat', 'tokens');
      fnum = str2num(regstrings{1}{1});
      fnums = [fnums fnum];
      dim = str2num(regstrings{1}{2});
      dims = [dims dim];
      evals = [evals str2num(regstrings{1}{3})];
      i = i+1;
    end
    tline = fgetl(fid);
  end
  fclose(fid);

  ufnums = unique(fnums);
  udims = unique(dims);

  load(files{1}, 'surrogateOpts');
  trainEvals = surrogateOpts.saveModelTrainingData;
  % trainEvals = [ 10 25 50 100 200 300 470 700 900 1200 1500 2000 2400 ];

  % 1st dim: functions
  % 2nd dim: trainEvals
  % 3rd dim: dimensions
  resRmse = nan(length(trainEvals), length(ufnums), length(udims));
  resCorr = nan(length(trainEvals), length(ufnums), length(udims));

  for i = 1:length(files)
    load(files{i}, 'BD', 'diagD', 'evalsReached', 'kendall', 'lambda', 'rmse', 'sigma', 'surrogateOpts', 'testsetX', 'testsetY', 'trainsetX', 'trainsetY', 'xmean');
    md = ModelFactory.createModel(modelStr, modelOpts, xmean');
    mdt = md.train(trainsetX, trainsetY, xmean', evals(i));
    yPredict = mdt.predict(testsetX);

    evalidx = find(trainEvals == evals(i));
    fnumidx = find(ufnums == fnums(i));
    dimidx = find(udims == dims(i));
    resRmse(evalidx, fnumidx, dimidx) = frmse(yPredict - testsetY);
    resCorr(evalidx, fnumidx, dimidx) = corr(yPredict, testsetY, 'type', 'kendall');
  end
end
