function MParams = InitModelParameters(ModelType,N,withModelOptimization)

if (withModelOptimization == 1)     MParams.MaxTrainingPoints = 2*(40 + floor(4 * (N^1.7)));
else                                MParams.MaxTrainingPoints = 40 + floor(4 * (N^1.7)) ;   end;

if (ModelType == 1) % RankSVM
    MParams.init = 0;

    % set the default hyper-parameters for surrogate-modelling here ...
    MParams.def_coeff(1) = 40 + floor(4 * (N^1.7)) ;      % number of training points
    MParams.def_coeff(2) = 6;                             % Cval = 10^coeff(2);
    MParams.def_coeff(3) = 3;                             % Ci(z) = Cval*((nAlpha-z).^coeff(3) );      
    MParams.def_coeff(4) = 1;                             % sigmaA
%        MParams.def_coeff(4) = 0;                             % sigmaA
        MParams.def_coeff(5) = 1000;
 %   MParams.def_coeff(5) = 0;                             % sigmaA
    MParams.coeff = MParams.def_coeff;
    
    % if you optimize the hyper-parameters of the model, please set Min Max ranges
    % and do not forget to update the BuildModel_RANKSVM.m  !
    MParams.xmin(1) = 4 * N;        MParams.xmax(1) = 2*(40 + floor(4 * (N^1.7)));   % number of training points
    MParams.xmin(2) = 0;            MParams.xmax(2) = 10;                            % Cval = 10^coeff(2);
    MParams.xmin(3) = 0;            MParams.xmax(3) = 6;                             % Ci(z) = Cval*((nAlpha-z).^coeff(3));  
    MParams.xmin(4) = 0.5;      	MParams.xmax(4) = 2;                             % sigmaA
        MParams.xmin(5) = 100;           MParams.xmax(5) = 1500;
%        MParams.xmin(4) = -10;           MParams.xmax(4) = 10;
  %  MParams.xmin(5) = -4;           MParams.xmax(5) = 4;  
    
    if (withModelOptimization == 1)     MParams.MaxTrainingPoints = MParams.xmax(1);
    else                                MParams.MaxTrainingPoints = MParams.def_coeff(1);   end;
    
end;

if (ModelType == 2) % SVR
    MParams.init = 0;

    % set the default hyper-parameters for surrogate-modelling here ...
    MParams.def_coeff(1) = 40 + floor(4 * (N^1.7)) ;      % number of training points
    MParams.def_coeff(2) = 6;                             % Cval = 10^coeff(2);
    MParams.def_coeff(3) = 3;                             % Ci(z) = Cval*((nAlpha-z).^coeff(3) );      
    MParams.def_coeff(4) = 1;                             % sigmaA
    
 %   MParams.def_coeff(5) = 1;                             % (1:nTrain).^coeff(5)
    
 %   MParams.def_coeff(6) = -5;                            % (ftrain(2) - ftrain(1))*( 10^coeff(6));
%    MParams.def_coeff(6) = 0;                             % sigmaA
    MParams.coeff = MParams.def_coeff;
    
    % if you optimize the hyper-parameters of the model, please set Min Max ranges
    % and do not forget to update the BuildModel_RANKSVM.m  !
    MParams.xmin(1) = 4 * N;        MParams.xmax(1) = 2*(40 + floor(4 * (N^1.7)));   % number of training points
    MParams.xmin(2) = 0;            MParams.xmax(2) = 10;                            % Cval = 10^coeff(2);
    MParams.xmin(3) = 0;            MParams.xmax(3) = 6;                             % Ci(z) = Cval*((nAlpha-z).^coeff(3));  
    MParams.xmin(4) = 0.5;      	MParams.xmax(4) = 2;                             % sigmaA
    
%    MParams.xmin(5) = 0.5;      	MParams.xmax(5) = 2;                             % (1:nTrain).^coeff(5)
%    MParams.xmin(6) = -5;           MParams.xmax(6) = 0;                            % (ftrain(2) - ftrain(1))*( 10^coeff(6));
%    MParams.xmin(6) = -2;           MParams.xmax(6) = 2;
    
    if (withModelOptimization == 1)     MParams.MaxTrainingPoints = MParams.xmax(1);
    else                                MParams.MaxTrainingPoints = MParams.def_coeff(1);   end;
    
end;


if (ModelType == 3) % RankStructSVM
    MParams.init = 0;

    % set the default hyper-parameters for surrogate-modelling here ...
    MParams.def_coeff(1) = 40 + floor(4 * (N^1.7)) ;      % number of training points
    MParams.def_coeff(2) = 6;                             % Cval = 10^coeff(2);
    MParams.def_coeff(3) = 3;                             % Ci(z) = Cval*((nAlpha-z).^coeff(3) );      
    MParams.def_coeff(4) = 0;                             % sigmaA
    MParams.coeff = MParams.def_coeff;
    
    % if you optimize the hyper-parameters of the model, please set Min Max ranges
    % and do not forget to update the BuildModel_RANKSVM.m  !
    MParams.xmin(1) = 4 * N;        MParams.xmax(1) = 2*(40 + floor(4 * (N^1.7)));   % number of training points
    MParams.xmin(2) = 0;            MParams.xmax(2) = 10;                            % Cval = 10^coeff(2);
    MParams.xmin(3) = 0;            MParams.xmax(3) = 6;                             % Ci(z) = Cval*((nAlpha-z).^coeff(3));  
    MParams.xmin(4) = 0;      	MParams.xmax(4) = 1;                             % sigmaA
    
    if (withModelOptimization == 1)     MParams.MaxTrainingPoints = MParams.xmax(1);
    else                                MParams.MaxTrainingPoints = MParams.def_coeff(1);   end;
    
end;