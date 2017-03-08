% script for creating matfile with algorithms to compare using variables
% generated in generatePlots_ppsn2016

funcSet.BBfunc = 1:24;
funcSet.dims = [2 3 5 10 20];

% CMA-ES
algorithm(1).name = 'CMA-ES';
algorithm(1).data = cmaesData;
algorithm(1).BBfunc = funcSet.BBfunc;
algorithm(1).dims = funcSet.dims;
algorithm(1).settings.PopSize = '(4 + floor(3*log(N)))';
algorithm(1).settings.Restarts = 4;
algorithm(1).color = getAlgColors('cmaes'); % [ 22,  22, 138]

% S-CMA-ES
algorithm(2).name = 'S-CMA-ES';
algorithm(2).data = genData;
algorithm(2).BBfunc = funcSet.BBfunc;
algorithm(2).dims = funcSet.dims;
algorithm(2).settings = gen_settings{genId};
algorithm(2).color = getAlgColors('scmaes'); % [178,  34,  34]

% DTS-CMA-ES
algorithm(3).name = 'DTS-CMA-ES';
algorithm(3).data = sd2Data_05_2pop;
algorithm(3).BBfunc = funcSet.BBfunc;
algorithm(3).dims = funcSet.dims;
algorithm(3).settings = sd2_r05_2pop_settings{sd2_r05_2pop_Id};
algorithm(3).color = getAlgColors('dtscmaes'); % [154, 205,  50]

% BIPOP-s*ACM-ES
algorithm(4).name = 'BIPOP-{}^{s*}ACMES-k';
algorithm(4).data = bbobDataReady(fullfile('exp', 'experiments', 'BIPOP-saACM-k'), funcSet);
algorithm(4).BBfunc = funcSet.BBfunc;
algorithm(4).dims = funcSet.dims;
algorithm(4).settings = 'Downloaded from http://coco.gforge.inria.fr/data-archive/2013/BIPOP-saACM-k_loshchilov_noiseless.tgz';
algorithm(4).color = getAlgColors('saacmes'); % [100, 149, 237]

% SMAC
algorithm(5).name = 'SMAC';
algorithm(5).data = bbobDataReady(fullfile('exp', 'experiments', 'SMAC'), funcSet);
algorithm(5).BBfunc = funcSet.BBfunc;
algorithm(5).dims = funcSet.dims;
algorithm(5).settings = 'Downloaded from http://coco.gforge.inria.fr/data-archive/2013/SMAC-BBOB_hutter_noiseless.tgz';
algorithm(5).color = getAlgColors('smac'); % [255, 155,   0]

% lmm-CMA-ES
algorithm(6).name = 'lmm-CMA-ES';
algorithm(6).data = bbobDataReady(fullfile('exp', 'experiments', 'lmm-CMA-ES'), funcSet);
algorithm(6).BBfunc = funcSet.BBfunc;
algorithm(6).dims = funcSet.dims;
algorithm(6).settings = 'Downloaded from http://coco.gforge.inria.fr/data-archive/2013/lmm-CMA-ES_auger_noiseless.tgz';
algorithm(6).color = getAlgColors('lmmcmaes'); % [255, 225,   0]

%% save algorithms data

resultname = which(mfilename);
save([resultname(1:end-1), 'mat'], 'algorithm');