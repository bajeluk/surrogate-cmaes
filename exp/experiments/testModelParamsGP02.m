function testModelParamsGP02(expPath, exp_id, varargin)
  param(1).name = 'covFcn';
  param(1).values = { 'covSEiso',  '{@covMaterniso, 3}', '{@covMaterniso, 5}' };
  param(1).fname  = { 'cov=SEiso', 'cov=Matern3'       , 'cov=Matern5'        };
  param(2).name = 'hyp.inf';
  param(2).values = { log(1e-3),  log(5e-3),  log(1e-2)  };
  param(2).fname  = { 'inf=1e-3', 'inf=5e-3', 'inf=1e-2' };
  param(3).name = 'hyp.lik';
  param(3).values = { log(0.01),  log(0.1),  log(1)  };
  param(3).fname  = { 'lik=0.01', 'lik=0.1', 'lik=1' };
  param(4).name = 'hyp.cov';
  param(4).values = { log([0.01; 0.05]), log([0.05; 0.25]), log([0.25; 1]), log([1; 5]) };
  param(4).fname  = { 'cov=0.01+0.05',   'cov=0.05+0.25',   'cov=0.25+1',   'cov=1+5'   };
  % THIS IS A HACK!!!
  SEisovalues =     { log([0.01; 0.1]), log([0.1; 10]), log([1; 1e3]), log([5; 1e4]) };
  SEisofnames =     { 'cov=0.01+0.1',   'cov=0.1+10',   'cov=1+1e3',   'cov=5+1e4'   };

  % file-path to the 'experiments' directory (here, new directory will
  % be placed and 'modeltrain_gene_02' directory is expected
  if (isempty(expPath))
    expPath = fileparts(mfilename('fullpath'));
  end

  nParams = length(param);
  for i = 1:nParams
      nParamValues(i) = length(param(i).values);
  end
  nCombinations = prod(nParamValues);

  % Combinations to try, can be specified as (a) parameter(s)
  switch nargin
  case 2
    combinationsToTry = 0:(nCombinations-1);
  case 3
    % combinations specified as vararg#1
    combinationsToTry = varargin{1};
  case 4
    % vararg#1 -- machine_id
    % vararg#2 -- total # of machines
    nMachines = varargin{2};
    assert(varargin{1} <= nMachines, '# of machines has to be > id');
    combsToMachine = floor(nCombinations / nMachines);
    startId = (varargin{1}-1) * combsToMachine;
    if (varargin{1} == nMachines)
      endId = nCombinations - 1;
    else
      endId = startId + combsToMachine - 1;
    end
    combinationsToTry = startId:endId;
    fprintf('Combinations to try: %d -- %d\n', startId, endId);
  end

  for i = combinationsToTry
      exactParamId = i; 
      filename = [expPath filesep exp_id filesep 'GP'];
      for j = 1:nParams
          ParamId = mod(exactParamId,nParamValues(j));
          % eval(['modelOpts.',param(j).name,' = ',param(j).values{ParamId+1},';'])
          fieldnames = strsplit(param(j).name, '.');
          switch length(fieldnames)
          case 1
            modelOpts.(param(j).name) = param(j).values{ParamId+1};
          case 2
            modelOpts.(fieldnames{1}).(fieldnames{2}) = param(j).values{ParamId+1};
          case 3
            modelOpts.(fieldnames{1}).(fieldnames{2}).(fieldnames{3}) = param(j).values{ParamId+1};
          otherwise
            % join the names with underscores
            modelOpts.(strjoin(fieldnames,'_')) = param(i).values{paramId+1};
          end
          fieldtag = param(j).fname{ParamId+1};

          % THIS IS A HACK!!!
          if (j == 4 && strcmp(modelOpts.covFcn, 'covSEiso'))
            fieldtag = SEisofnames{ParamId+1};
            modelOpts.hyp.cov = SEisovalues{ParamId+1};
          end

          filename = [filename '_' param(j).fname{ParamId+1}];
          exactParamId = (exactParamId-ParamId)/nParamValues(j);
      end
      disp(modelOpts)
      if (isfield(modelOpts, 'hyp'))
        disp(modelOpts.hyp);
        disp(modelOpts.hyp.cov);
      end

      % 
      
      trainlist = [expPath filesep 'modeltrain_gene_02' filesep 'modeltrain_trainlist.txt'];
      [RMSE,kor,resTime,resLiks,resErrs] = modelTrainTest('gp',modelOpts,trainlist);
      
      [~,mess,messid] = mkdir([expPath filesep exp_id]);
      save([filename '.mat'], 'RMSE', 'kor', 'resTime', 'resLiks', 'resErrs', 'filename'); 
  end
end
