paramFolder = fullfile('exp','experiments','modeltrain_rf_01'); %'RFparams');

fileList = dir([paramFolder filesep '*.mat']);

fileIdxs = 1:length(fileList);
% chosen good GP model settings:
% fileIdxs = [3 4 7 8 11 23 32 39 40 43 44 47 48 63 75 76 79 87 88 95 99]; 

%% Loading files

Dim2= [];Dim5= [];Dim10= [];
line = 1;
for i = fileIdxs
    fprintf('%d  %s\n',i,fileList(i).name);
    loadKor = load([paramFolder filesep fileList(i).name],'kor');
    kor = loadKor.kor;
    Dim2(:,:,line) = loadKor.kor(:,:,1);
    Dim5(:,:,line) = loadKor.kor(:,:,2);
    Dim10(:,:,line) = loadKor.kor(:,:,3);
    line = line + 1;
end

%% Drawing

% draw boxplots and bars of number of non NaN values
numOfFunc = size(Dim2,2);
funcLabels = {'f1','f2','f3','f5','f6','f8','f10','f11','f12','f13','f14','f20','f21'};
for i = 1:numOfFunc
    figure('Name',funcLabels{i})
    
    rDim2 = reshape(Dim2(:,i,:),size(Dim2,1),size(Dim2,3));
    subplot(6,1,1);
    boxplot(rDim2);
    subplot(6,1,2);
    bar(sum(~isnan(rDim2),1));
    
    rDim5 = reshape(Dim5(:,i,:),size(Dim2,1),size(Dim2,3));
    subplot(6,1,3);
    boxplot(rDim5);
    subplot(6,1,4);
    bar(sum(~isnan(rDim5),1));    
    
    rDim10 = reshape(Dim10(:,i,:),size(Dim2,1),size(Dim2,3));
    subplot(6,1,5);
    boxplot(rDim10);
    subplot(6,1,6);
    bar(sum(~isnan(rDim10),1));
end
% 
% figure('Name','f2')
% subplot(3,1,1);
% boxplot(kor2(:,1:36));
% subplot(3,1,2);
% boxplot(kor2(:,37:72));
% subplot(3,1,3);
% boxplot(kor2(:,73:108));
% 
% figure('Name','f8')
% subplot(3,1,1);
% boxplot(kor3(:,1:36));
% subplot(3,1,2);
% boxplot(kor3(:,37:72));
% subplot(3,1,3);
% boxplot(kor3(:,73:108));
% 
% figure('Name','f10')
% subplot(3,1,1);
% boxplot(kor4(:,1:36));
% subplot(3,1,2);
% boxplot(kor4(:,37:72));
% subplot(3,1,3);
% boxplot(kor4(:,73:108));

clear paramFolder fileList i loadKor
