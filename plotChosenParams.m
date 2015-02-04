% plotChosenParams
close all

dims = [2 5 10];
funcLabels = {'f1','f2','f3','f5','f6','f8','f10','f11','f12','f13','f14','f20','f21'};
functionLabels = {'Sphere','Ellipsoid separable','Rastrigin','Linear slope','Attractive sector function',...
                  'Rosenbrock','Ellipsoid with monotone x-transformation',...
                  'Discus with monotone x-transformation','Bent cigar',...
                  'Sharp ridge','Sum of different powers','Schwefel','Gallagher'};
models = {'gp','rf'};
 
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

modelLabels={'GPb','GPm','RFb','RFm'};

% draw boxplots of RMSE and correlation for every function and dimension
for f = 1:nFunc
    
    % plot correlation
    figure('Name',['f',num2str(f),' corr'],'Units','centimeters','Position',[f/10,15-f,15 5]);  
    for D = 1:nDim
          subplot(1, nDim, D); % (f-1)*nDim + D
          boxplot([gpkor(:,f,D,bestCorr_gp(f,D)),gpkor(:,f,D,medCorr_gp(f,D)),...
              rfkor(:,f,D,bestCorr_rf(f,D)),rfkor(:,f,D,medCorr_rf(f,D))],modelLabels);
          title([int2str(dims(D)),' D']);
          ylim([-0.5 1.05]);
          if D == 1
               ylabel('Correlation');
          end
          set(gca,'Position',[0.1+0.3*(D-1) 0.15 0.25 0.7])
    end
    
    % plot RMSE
    figure('Name',['f',num2str(f),' RMSE'],'Units','centimeters','Position',[17+f/10,15-f,15,5]);
    for D = 1:nDim
          subplot(1, nDim, D); % (f-1)*nDim + D
          boxplot(log([gpRMSE(:,f,D,bestRMSE_gp(f,D)),gpRMSE(:,f,D,medRMSE_gp(f,D)),...
              rfRMSE(:,f,D,bestRMSE_rf(f,D)),rfRMSE(:,f,D,medRMSE_rf(f,D))]),modelLabels);
          title([int2str(dims(D)),'D']);
          ylim([-20 2]);
          if D == 1
              ylabel('log RMSE');
          end
          set(gca,'Position',[0.1+0.3*(D-1) 0.15 0.25 0.7])
    end
end