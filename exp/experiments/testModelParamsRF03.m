function testModelParamsRF03(path, exp_id)
  param(1).name = 'nTrees';
  param(1).values = {100};
  param(2).name = 'nBestPoints';
  param(2).values = {5};
  param(3).name = 'minLeaf';
  param(3).values = {2,5,8};
  param(4).name = 'inputFraction';
  param(4).values = {1};

  nParams = length(param);
  for i = 1:nParams
      nParamValues(i) = length(param(i).values);
  end
  nCombinations = prod(nParamValues);

  for i = 0:nCombinations - 1
      exactParamId = i; 
      filename = [path filesep exp_id filesep 'RF'];
      for j = 1:nParams
          ParamId = mod(exactParamId,nParamValues(j));
          eval(['modelOpts.',param(j).name,' = ',num2str(param(j).values{ParamId+1}),';'])
          if strcmp(param(j).name,'inputFraction')
              filename = [filename,'_',param(j).name,'_',num2str(10*param(j).values{ParamId+1})];
          else
              filename = [filename,'_',param(j).name,'_',num2str(param(j).values{ParamId+1})];
          end
          exactParamId = (exactParamId-ParamId)/nParamValues(j);
      end
      disp(modelOpts)
      
      trainlist = [path filesep 'modeltrain_gene_02' filesep 'modeltrain_trainlist.txt'];
      [RMSE,kor,resTime] = modelTrainTest('rf',modelOpts,trainlist);
      
      mkdir([path filesep exp_id]);
      save([filename,'.mat'],'RMSE','kor','resTime'); 
  end
end
