% plotChosenParams

models = {'gp','rf'};
dims = [2 5 10];
 
% gather all data
for mod = models
     m = mod{1};
     Folder = fullfile('exp','experiments',['modeltrain_',m,'_01']);
     fileList = dir([Folder filesep '*.mat']);
     nFiles = length(fileList);
     % process each file
     for i = 1:nFiles 
%             fprintf('%d  %s\n',i,fileList(i).name);
            Data = load([Folder filesep fileList(i).name],'RMSE','kor');
            RMSE(:,:,:,i) = Data.RMSE;
            kor(:,:,:,i) = Data.kor;
     end
     eval([m,'RMSE = RMSE;']);
     eval([m,'kor = kor;']);
     
     % find best and median parametres
     [bestRMSE,medRMSE,bestCorr,medCorr] = findBestAndAverageParams(Folder);
     
     eval(['bestRMSE_',m,' = bestRMSE;']);
     eval(['medRMSE_',m,' = medRMSE;']);
     eval(['bestCorr_',m,' = bestCorr;']);
     eval(['medCorr_',m,' = medCorr;']);
end

% matrix coordinates: 1 tresholds, 2 functions, 3 dimensions, 4 parameter settings
[~, nFunc, nDim, ~] = size(gpRMSE); 

% draw boxplots of RMSE and correlation for every function and dimension
for f = 1:1 %nFunc
    
    % plot correlation
    figure();
    for D = 1:nDim
          subplot(1, nDim, D); % (f-1)*nDim + D
          boxplot([gpkor(:,f,D,bestCorr_gp(f,D)),gpkor(:,f,D,medCorr_gp(f,D)),rfkor(:,f,D,bestCorr_rf(f,D)),rfkor(:,f,D,medCorr_rf(f,D))]);
          title(['f',int2str(f),' ',int2str(dims(D)),'D']);
          ylim([-0.5 1.05]);
    end
    
    % plot RMSE
    figure();
    for D = 1:nDim
          subplot(1, nDim, D); % (f-1)*nDim + D
          boxplot(log([gpRMSE(:,f,D,bestRMSE_gp(f,D)),gpRMSE(:,f,D,medRMSE_gp(f,D)),rfRMSE(:,f,D,bestRMSE_rf(f,D)),rfRMSE(:,f,D,medRMSE_rf(f,D))]));
          title(['f',int2str(f),' ',int2str(dims(D)),'D']);
%           ylim([-0.5 1.05]);
    end
end