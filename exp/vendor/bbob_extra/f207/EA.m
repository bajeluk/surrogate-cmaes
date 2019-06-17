function [  ] = EA()
%EA Summary of this function goes here
%   Detailed explanation goes here


%% Constant setting
numFrequencies = 1; % number of frequencies
dominant_frequency = 0.7; % dominant frequency, used for single frequency
numSphere = 5; % number of buoys

width = 100; % width of farm
height = 100;  % height of farm
radiusSphere = 2*ones(1, numSphere);

%% print head
test_name = int2str(feature('getpid'));
fprintf('PID: %s.\n', test_name);
fprintf('Number of Buoys: %d.\n', numSphere);
fprintf('Number of Frequencies: %d.\n', numFrequencies);

%% initialising
% loop until the feasiable instance
rng('shuffle');
notdone = true;
while notdone
    xSphere = randi(width, 1, numSphere);
    ySphere = randi(height, 1, numSphere);
    
    notdone = ~layoutConstraint(xSphere, ySphere);
end


%% save the instance
instance_name = ['000-' datestr(datetime(),'yyyymmddTHHMMSS')];
folder_name = strcat('layouts/mutation_', test_name,'/');
[~,~,~] = mkdir(folder_name);
save(strcat(folder_name, instance_name, '.mat'), 'xSphere', 'ySphere');

% load('layouts/mutation/301-20160108T151520.mat');
% instance_name = '301-20160108T151520';

fprintf('G#Instance\tbestRCW\tbestPower\tRCW\tPower\tTime(s)\n');

tic;
if numFrequencies == 1
    warning('off','all');
    [MlRCWWeighted, MlParrayAverage, ParrayBuoy] = arrayBuoyPlacementSingleFreq(radiusSphere, xSphere, ySphere, dominant_frequency, 0);
else
    gcp;
    pctRunOnAll warning('off','all');
    [MlRCWWeighted, MlParrayAverage] = arrayBuoyPlacement_v2(radiusSphere, xSphere, ySphere, numFrequencies);
end
    
MlTime = toc;


% scatter3(xSphere,ySphere,ParrayBuoy,100,ParrayBuoy,'filled');
% view(2);
% colormap('jet');
% colorbar;



bestMltRCW = MlRCWWeighted;
bestMltPower = MlParrayAverage;


fprintf('%s\t%.4g\t%.5g\t%.4g\t%.5g\t%.2f\n',instance_name,bestMltRCW,bestMltPower,MlRCWWeighted,MlParrayAverage,MlTime);


%% iterating
for i=1:400

    % mutation
    selected = randi(numSphere);
    % store the current instance
    cx = xSphere(selected);
    cy = ySphere(selected);
    
    % loop until the feasiable mutation
    notdone = true;
    while notdone
        xSphere(selected) = randi(width);
        ySphere(selected) = randi(height);
        
        notdone = ~layoutConstraint(xSphere, ySphere);
    end
    
    % save the current instance
    instance_name = [sprintf('%03d',i) '-' datestr(datetime(),'yyyymmddTHHMMSS')];
    save(strcat(folder_name, instance_name, '.mat'), 'xSphere', 'ySphere');
    
%     figure;
%     scatter(xSphere,ySphere);
    
    tic;
    if numFrequencies == 1
        [MlRCWWeighted, MlParrayAverage] = arrayBuoyPlacementSingleFreq(radiusSphere, xSphere, ySphere, dominant_frequency, 0);
    else
        [MlRCWWeighted, MlParrayAverage] = arrayBuoyPlacement_v2(radiusSphere, xSphere, ySphere, numFrequencies);
    end
    MlTime = toc;

    %selection 
    if MlParrayAverage > bestMltPower
        bestMltRCW = MlRCWWeighted;
        bestMltPower = MlParrayAverage;
    else % return to current instance
        xSphere(selected) = cx;
        ySphere(selected) = cy;
    end
    
    fprintf('%s\t%.4g\t%.5g\t%.4g\t%.5g\t%.2f\n',instance_name,bestMltRCW,bestMltPower,MlRCWWeighted,MlParrayAverage,MlTime);
    
end

% save the best instance
save(strcat(folder_name, 'best-',datestr(datetime(),'yyyymmddTHHMMSS'),'.mat'), 'xSphere', 'ySphere');
% delete(gcp);


end

