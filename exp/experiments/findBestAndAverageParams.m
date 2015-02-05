function [bestRMSE,medianRMSE,bestCorr,medianCorr] = findBestAndAverageParams(dataFolder)

%     RFFolder = fullfile('exp','experiments','modeltrain_rf_01'); % RF params
%     GPFolder = fullfile('exp','experiments','modeltrain_gp_01'); % GP params

    % load file-names in the directory dataFolder
    fileList = dir([dataFolder filesep '*.mat']);
    
    nFiles = length(fileList);
    % process each file
    for i = 1:nFiles 
%         fprintf('%d  %s\n',i,fileList(i).name);
        Data = load([dataFolder filesep fileList(i).name],'RMSE','kor');
        RMSE = Data.RMSE;
        kor = Data.kor;
        
        dims = size(RMSE,3);
        
        for dim = 1:dims
             % identify non-NaN entries
             validRMSE = ~isnan(RMSE(:,:,dim));
             validCorr = ~isnan(kor(:,:,dim));
             % calculate mean of correlation and RMSE for each function through eval. thresholds
             for j = 1:size(validRMSE,2)
                 medRMSE(j,dim,i) = median(RMSE(validRMSE(:,j),j,dim));
                 medCorr(j,dim,i) = median(kor(validCorr(:,j),j,dim));
             end
        end
    end
    
    % find best correlation and RMSE
    [~,bestRMSE] =  min(medRMSE,[],3);
    [~,bestCorr] =  max(medCorr,[],3);
    
    % find median correlation and RMSE
    [~,medianCorr] = min(abs(repmat(median(medCorr,3),[1,1,nFiles])-medCorr),[],3);
    [~,medianRMSE] = min(abs(repmat(median(medRMSE,3),[1,1,nFiles])-medRMSE),[],3);
end

