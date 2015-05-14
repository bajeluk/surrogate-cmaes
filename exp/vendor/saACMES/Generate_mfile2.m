
for iexp=1:60
    BIPOP = 1;
    newRestartRules = 0;
    noisy = 0;
    CMAactive = 1;
    modelType = 1;
    lambdaMult = 1;
    muMult = 1;
    
   % largeLambdaMinIter = 1;
    
    if ((iexp > 0) && (iexp <= 15))
        ind = iexp;
        part = 1;
  %      largeLambdaMinIter = 1;
    end;
    if ((iexp > 15) && (iexp <= 30))
        ind = iexp - 15;
        part = 2;
 %       largeLambdaMinIter = 2;
    end;
    if ((iexp > 30) && (iexp <= 45))
        ind = iexp - 30;
        part = 3;
%        largeLambdaMinIter = 3;
    end;
    if ((iexp > 45) && (iexp <= 60))
        ind = iexp - 45;
        part = 4;
%        largeLambdaMinIter = 3;
    end;
    
    largeLambdaMinIter = 3;

    file = fopen(['exp' num2str(iexp) '.m'],'w');
    
    inststr = num2str(ind);
    
    fprintf(file, ['global settings;\n']);
    fprintf(file, ['\n']);
    fprintf(file, ['%problem A\n']);
    
    

    if (part == 1)     
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [20];\n']);                  
        fprintf(file, ['settings.funs = [1,2,5:15];\n']);      
    end;
    if (part == 2)   
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [20];\n']);                  
        fprintf(file, ['settings.funs = [16,17,18,19];\n']);  
    end;
    if (part == 3)   
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [20];\n']);                  
        fprintf(file, ['settings.funs = [20,21,22];\n']);  ;
    end;
    if (part == 4)   
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [20];\n']);                  
        fprintf(file, ['settings.funs = [23,24,3,4];\n']);
    end;


    fprintf(file, ['settings.pathname = ''d' num2str(iexp) ''';\n']);
    fprintf(file, ['settings.algname = '''';\n']);
    fprintf(file, ['settings.ntarray = [1];\n']);
    fprintf(file, ['settings.savfile = ''r' num2str(iexp) ''';\n']);
    fprintf(file, ['\n']);
    
    
    fprintf(file, ['settings.BIPOP = ' num2str(BIPOP) '; \n']);
    fprintf(file, ['settings.newRestartRules = ' num2str(newRestartRules) '; \n']);
    fprintf(file, ['settings.noisy = ' num2str(noisy) ';\n']);
    fprintf(file, ['settings.CMAactive = ' num2str(CMAactive) ';\n']);
    fprintf(file, ['settings.withFileDisp = 1;\n']);
    fprintf(file, ['settings.withSurr = 1;\n']);
    fprintf(file, ['settings.modelType = ' num2str(modelType) ';\n']);
    fprintf(file, ['settings.withModelEnsembles = 0;\n']);
    fprintf(file, ['settings.withModelOptimization = 1;\n']);
    fprintf(file, ['settings.hyper_lambda = 20;\n']);
    fprintf(file, ['settings.iSTEPminForHyperOptimization = 1;\n']);
    fprintf(file, ['settings.MaxEvals = ''1e6*dim'';\n']);
    fprintf(file, ['settings.MaxEvalsWithSurrogate = ''1e4*20'';\n']);
    fprintf(file, ['settings.lambdaMult = ' num2str(lambdaMult) ';\n']);
    fprintf(file, ['settings.muMult = ' num2str(muMult) ';\n']);
    fprintf(file, ['settings.largeLambdaMinIter = ' num2str(largeLambdaMinIter) ';\n']);

    fprintf(file, ['\n']);
    fprintf(file, ['settings.withDisp = 0;\n']);
    fprintf(file, ['settings.maxStepts = 20;\n']);
    fprintf(file, ['settings.maxerr = 0.45;\n']);
    fprintf(file, ['settings.alpha = 0.20;\n']);
    fprintf(file, ['settings.iterstart = 10;\n']);

    fprintf(file, ['Adapter();\n']);

    
    fclose(file);
    
end;

z = 0;