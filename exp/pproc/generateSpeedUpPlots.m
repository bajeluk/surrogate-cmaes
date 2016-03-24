function generateSpeedUpPlots()
  % function for making speed up graphs on gecco 2015 abstract
  % The speed up is considered between the CMA-ES and GP, or RF
  
  close all

  exppath = fullfile('exp', 'experiments');
  gppath = fullfile(exppath, 'exp_geneEC_06_all', 'exp_geneEC_06');
  rfpath = {[gppath,'_rflite'], [gppath,'_rf5_2'], [gppath,'_rf5_3'], [gppath,'_rf5_4']};
  rfpath = {[gppath,'_rf5_2'], [gppath,'_rf5_3'], [gppath,'_rf5_4']};
  cmaespath = fullfile(gppath, 'cmaes_results');
  plotResultsFolder = '/tmp';

  funcSet.BBfunc = [1,2,3,5,6,8,10,11,12,13,14,20,21];
  funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
  funcSet.dims = [2,5,10];
  funcSet.dimsInv = inverseIndex(funcSet.dims);
  
  rf_evals = dataReady(rfpath, funcSet, 1, 'rflite');
  gp_evals = dataReady(gppath, funcSet, 6);
  cmaes_evals = dataReady(cmaespath, funcSet, 1, 'cmaes');
  
  hRFall = speedUpPlot(rf_evals,cmaes_evals,'cmaes',funcSet);
  pdfname = fullfile('speedUpRF');
  print2pdf(hRFall,pdfname,1)

%   speedUpPlot(gp_evals,cmaes_evals,'cmaes',funcSet);
%   
  funcSet.BBfunc = [1,2,3,5,6,8];
  h1 = speedUpPlotCompare(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);
  
  funcSet.BBfunc = [10,11,12,14,20,21];
  h2 = speedUpPlotCompare(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);
  
  pdfname = fullfile(plotResultsFolder,'speedUp');
  print2pdf([h1 h2],{[pdfname,'A.pdf'],[pdfname,'B.pdf']},1)

end