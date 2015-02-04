function [resRmse, resCorr, resTime, varargout] = modelTrainTest(modelStr, modelOpts, filelist)
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
  pathList = fileparts(filelist);
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
      % collect function numbers, dimensions and #fevals
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

  load([pathList filesep files{1}], 'surrogateOpts');
  trainEvals = surrogateOpts.saveModelTrainingData;
  % trainEvals = [ 10 25 50 100 200 300 470 700 900 1200 1500 2000 2400 ];

  % Output arrays:
  %   1st dim: functions
  %   2nd dim: trainEvals
  %   3rd dim: dimensions
  resRmse = nan(length(trainEvals), length(ufnums), length(udims));
  resCorr = nan(length(trainEvals), length(ufnums), length(udims));
  resTime = nan(length(trainEvals), length(ufnums), length(udims));
  if (strcmpi(modelStr, 'gp'))
    resLiks = nan(length(trainEvals), length(ufnums), length(udims));
    resErrs = nan(length(trainEvals), length(ufnums), length(udims));
  end

  for i = 1:length(files)
    load([pathList filesep files{i}], 'BD', 'diagD', 'evalsReached', 'kendall', 'lambda', 'rmse', 'sigma', 'surrogateOpts', 'testsetX', 'testsetY', 'trainsetX', 'trainsetY', 'xmean');
    disp(['.. testing file ' files{i} ' ..']);
    md = ModelFactory.createModel(modelStr, modelOpts, xmean');
    isTrained = false;
    j = 0;
    t = tic();
    while (~isTrained && j < 5)
      mdt = md.train(trainsetX, trainsetY, xmean', evals(i));
      isTrained = mdt.isTrained();
      j = j+1;
    end
    timeSpent = toc(t);
    % indices of #evals, fnum and dimensionality of the test-case
    %
    % TODO: this does not count with multiple test-cases with the
    %   same f-num and dim!!! Rewrite this to calculate median/or mean!
    %
    evalidx = find(trainEvals == evals(i));
    fnumidx = find(ufnums == fnums(i));
    dimidx = find(udims == dims(i));
    if (strcmpi(modelStr, 'gp'))
      % record final the likelihood achieved by optimization
      resLiks(evalidx, fnumidx, dimidx) = mdt.trainLikelihood;
      resErrs(evalidx, fnumidx, dimidx) = j*mdt.nErrors;
    end

    if (isTrained)
      % predict on the test-set
      yPredict = mdt.predict(testsetX);
      % record the RMSE, Kendall's correlation and CPU time
      resRmse(evalidx, fnumidx, dimidx) = frmse(yPredict - testsetY);
      resCorr(evalidx, fnumidx, dimidx) = corr(yPredict, testsetY, 'type', 'kendall');
      resTime(evalidx, fnumidx, dimidx) = timeSpent;
      if (strcmpi(modelStr, 'gp'))
        quality = mdt.trainLikelihood;
      else
        quality = 0;
      end
      fprintf('.. rmse = %f ,  time = %0.2f s,  lik = %0.2f,  corr = %0.3f\n', resRmse(evalidx, fnumidx, dimidx), timeSpent, quality, resCorr(evalidx, fnumidx, dimidx));
    end
    % if (strcmpi(modelStr, 'gp'))
    %   fprintf('.. final hyp = ');
    %   disp(unwrap(mdt.hyp)');
    % end
  end

  if (strcmpi(modelStr, 'gp') && nargout > 3)
    % record final the likelihood achieved by optimization
    varargout(1) = {resLiks};
    varargout(2) = {resErrs};
  end
end
