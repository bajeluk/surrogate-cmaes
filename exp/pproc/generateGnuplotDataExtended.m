function generateGnuplotDataExtended(gnuplotFile, exp_results, exp_cmaes_results, maxfunevals)
  fid = fopen(gnuplotFile, 'w');
  [mBestf, q1Bestf, q3Bestf] = statisticsFromYEvals(exp_results.y_evals, maxfunevals, 1);
  [mRMSE, q1RMSE, q3RMSE] = statisticsFromYEvals(exp_results.y_evals, maxfunevals, 3);
  [mCorr, q1Corr, q3Corr] = statisticsFromYEvals(exp_results.y_evals, maxfunevals, 4);
  names = {'sigma', 'diagRatio', 'minstd', 'maxstd', 'testRangeErr', 'popRangeErr', 'nTrainPoints', 'reevalPointsDist'};
  stats = [];
  colNames = '';
  nCols = 0
  while (size(exp_results.y_evals{1}, 2) >= (5 + nCols))
    nCols = nCols + 1;
    stats(:,nCols) = statisticsFromYEvals(exp_results.y_evals, maxfunevals, 4+nCols);
    colNames = [colNames names{nCols} ' '];
  end

  [mCmaes, q1Cmaes, q3Cmaes] = statisticsFromYEvals(exp_cmaes_results.y_evals, maxfunevals, 1);

  fprintf(fid, '# evals DataMedian DataQ1 DataQ3 CmaesMedian CmaesQ1 CmaesQ3 ModelRMSE ModelCorr %s\n', colNames);
  for i = 1:maxfunevals
    str = sprintf('%d %e %e %e %e %e %e %e %e', i, mBestf(i), q1Bestf(i), q3Bestf(i), mCmaes(i), q1Cmaes(i), q3Cmaes(i), mRMSE(i), mCorr(i));
    for j = 1:nCols
      str = [str ' ' num2str(stats(i,j))];
    end
    fprintf(fid, '%s\n', str);
  end

  fclose(fid);
end

