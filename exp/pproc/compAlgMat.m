% script for creating matfile with algorithms to compare using variables
% generated in generatePlots_ppsn2016

dims = [2 3 5 10 20];

% CMA-ES
algorithm(1).name = 'CMA-ES';
algorithm(1).data = cmaesData;
algorithm(1).dims = dims;
algorithm(1).settings.PopSize = '(4 + floor(3*log(N)))';
algorithm(1).settings.Restarts = 4;
algorithm(1).color = [ 22,  22, 138];

% S-CMA-ES
algorithm(2).name = 'S-CMA-ES';
algorithm(2).data = genData;
algorithm(2).dims = dims;
algorithm(2).settings = gen_settings{genId};
algorithm(2).color = [178,  34,  34];

% DTS-CMA-ES
algorithm(3).name = 'DTS-CMA-ES';
algorithm(3).data = sd2Data_05_2pop;
algorithm(3).dims = dims;
algorithm(3).settings = sd2_r05_2pop_settings{sd2_r05_2pop_Id};
algorithm(3).color = [154, 205,  50];

% BIPOP-s*ACM-ES
algorithm(4).name = 'BIPOP-{}^{s*}ACMES-k';
algorithm(4).data = saacmesData;
algorithm(4).dims = dims;
algorithm(4).settings = 'Downloaded from http://coco.gforge.inria.fr/data-archive/2013/BIPOP-saACM-k_loshchilov_noiseless.tgz';
algorithm(4).color = [100, 149, 237];

% SMAC
algorithm(5).name = 'SMAC';
algorithm(5).data = smacData;
algorithm(5).dims = dims;
algorithm(5).settings = 'Downloaded from http://coco.gforge.inria.fr/data-archive/2013/SMAC-BBOB_hutter_noiseless.tgz';
algorithm(5).color = [255, 155,   0];

%% save algorithms data

resultname = which(mfilename);
save([resultname(1:end-1), 'mat'], 'algorithm');