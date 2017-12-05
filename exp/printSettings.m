function printSettings(fid, exp_settings, exp_results, surrogateParams, cmaesParams)
  fprintf(fid, '===== Experiment: %s =====\n\n', exp_settings.exp_id);
  fprintf(fid, '== BBOB experiment settings: ==\n');
  fprintf(fid, sprintfStruct(exp_settings));
  fprintf(fid, '%15s: %f\n', 'time elapsed', exp_results.time);
  fprintf(fid, '\n== Surrogate model parameters: ==\n');
  fprintf(fid, sprintfStruct(surrogateParams));
  fprintf(fid, '\n== CMA-ES parameters: ==\n');
  fprintf(fid, sprintfStruct(cmaesParams));
    fprintf(fid, '\n== CMA-ES surrogate model options: ==\n');
  fprintf(fid, sprintfStruct(surrogateParams.modelOpts));
  fprintf(fid, '\n== Numerical results: ==\n\n');
  fprintf(fid, 'fbests:\n%s\n\n', num2str(exp_results.fbests));
  fprintf(fid, 'f075:\n%s\n\n', num2str(exp_results.f075));
  fprintf(fid, 'f050:\n%s\n\n', num2str(exp_results.f050));
  fprintf(fid, 'f025:\n%s\n\n', num2str(exp_results.f025));
end