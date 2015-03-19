param(1).name = 'nTrees';
param(1).values = {100};% {100,400,1000};
param(2).name = 'nBestPoints';
param(2).values = {3,5}; % {0,1,3,5};
param(3).name = 'minLeaf';
param(3).values = {2,5,8}; % {2,5,8};
param(4).name = 'inputFraction';
param(4).values = {0.8,1};% {0.5,0.8,1};
param(5).name = 'transformCoordinates';
param(5).values = {true}; % {true,false};

nParams = length(param);
for i = 1:nParams
    nParamValues(i) = length(param(i).values);
end
nCombinations = prod(nParamValues);

for i = 0:nCombinations - 1
    exactParamId = i; 
    filename = 'exp/experiments/RFparams/RF';
    for j = 1:nParams
        ParamId = mod(exactParamId,nParamValues(j));
        eval(['modelOpts.',param(j).name,' = ',num2str(param(j).values{ParamId+1}),';'])
        switch param(j).name
            case 'inputFraction'
                filename = [filename,'_',param(j).name,'_',num2str(10*param(j).values{ParamId+1})];
            case 'transformCoordinates'
                filename = [filename,'_transCoor_',num2str(param(j).values{ParamId+1})];
            otherwise
                filename = [filename,'_',param(j).name,'_',num2str(param(j).values{ParamId+1})];
        end
        exactParamId = (exactParamId-ParamId)/nParamValues(j);
    end
    disp(modelOpts)
    
    trainlist = 'exp/experiments/modeltrain_gene_01/modeltrain_trainlist.txt';
    [RMSE,kor] = modelTrainTest('rf',modelOpts,trainlist);
    
    mkdir('exp/experiments/RFparams');
    save([filename,'.mat'],'RMSE','kor');
end

clear
