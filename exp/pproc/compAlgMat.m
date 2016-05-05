% script for creating matfile with algorithms to compare using variables
% generated in generatePlots_ppsn2016

dims = [2 3 5 10 20];

% CMA-ES
algorithm(1).name = 'cmaes';
algorithm(1).data = cmaesData;
algorithm(1).dims = dims;
algorithm(1).settings.PopSize = '(4 + floor(3*log(N)))';
algorithm(1).settings.Restarts = 4;

% S-CMA-ES
algorithm(2).name = 'scmaes';
algorithm(2).data = genData;
algorithm(2).dims = dims;
algorithm(2).settings = gen_settings{genId};

% DTS-CMA-ES
algorithm(3).name = 'dtscmaes';
algorithm(3).data = sd2Data_05_2pop;
algorithm(3).dims = dims;
algorithm(3).settings = sd2_r05_2pop_settings{sd2_r05_2pop_Id};

% BIPOP-s*ACM-ES
algorithm(4).name = 'saacmes';
algorithm(4).data = saacmesData;
algorithm(4).dims = dims;
algorithm(4).settings = 'Downloaded from http://coco.gforge.inria.fr/data-archive/2013/BIPOP-saACM-k_loshchilov_noiseless.tgz';

% SMAC
algorithm(5).name = 'smac';
algorithm(5).data = smacData;
algorithm(5).dims = dims;
algorithm(5).settings = 'Downloaded from http://coco.gforge.inria.fr/data-archive/2013/SMAC-BBOB_hutter_noiseless.tgz';

%% save algorithms data

resultname = which(mfilename);
save([resultname(1:end-1), 'mat'], 'algorithm');