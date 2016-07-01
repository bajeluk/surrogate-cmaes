function generateGnuplotData(gnuplotFile, exp_results, exp_cmaes_results, maxfunevals)
  fid = fopen(gnuplotFile, 'w');
  [mData, q1Data, q3Data] = statisticsFromYEvals(exp_results.y_evals, maxfunevals, 1);
  [mCmaes, q1Cmaes, q3Cmaes] = statisticsFromYEvals(exp_cmaes_results.y_evals, maxfunevals, 1);

  fprintf(fid, '# evals DataMedian DataQ1 DataQ3 CmaesMedian CmaesQ1 CmaesQ3\n');
  for i = 1:maxfunevals
    fprintf(fid, '%d %e %e %e %e %e %e\n', i, mData(i), q1Data(i), q3Data(i), mCmaes(i), q1Cmaes(i), q3Cmaes(i));
  end

  fclose(fid);
end

