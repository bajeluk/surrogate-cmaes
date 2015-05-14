
for iexp=1:165
    BIPOP = 1;
    newRestartRules = 0;
    noisy = 0;
    CMAactive = 1;
    modelType = 1;
    lambdaMult = 1;
    muMult = 1;
    
    if ((iexp > 0) && (iexp <= 15))
        ind = iexp;
        part = 1;
    end;
    if ((iexp > 15) && (iexp <= 30))
        ind = iexp - 15;
        part = 2;
    end;
    if ((iexp > 30) && (iexp <= 45))
        ind = iexp - 30;
        part = 3;
    end;
    if ((iexp > 45) && (iexp <= 60))
        ind = iexp - 45;
        part = 4;
    end;
    if ((iexp > 60) && (iexp <= 75))
        ind = iexp - 60;
        part = 5;
    end;
    if ((iexp > 75) && (iexp <= 90))
        ind = iexp - 75;
        part = 6;
    end;
    if ((iexp > 90) && (iexp <= 105))
        ind = iexp - 90;
        part = 7;
    end;
    if ((iexp > 105) && (iexp <= 120))
        ind = iexp - 105;
        part = 8;
    end;
    if ((iexp > 120) && (iexp <= 135))
        ind = iexp - 120;
        part = 9;
    end;
    if ((iexp > 135) && (iexp <= 150))
        ind = iexp - 135;
        part = 10;
    end;
    if ((iexp > 150) && (iexp <= 165))
        ind = iexp - 150;
        part = 11;
    end;


    
    file = fopen(['exp' num2str(iexp) '.m'],'w');
    
    inststr = num2str(ind);
    
    fprintf(file, ['global settings;\n']);
    fprintf(file, ['\n']);
    fprintf(file, ['%problem A\n']);
    
    
    if (part == 1)     
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [10,20];\n']);                  
        fprintf(file, ['settings.funs = [1,2,5:19];\n']);      
    end;
    if (part == 2)          
        if (ind == 1)   fprintf(file, ['settings.dims = [10];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [1];\n']);     end;
        if (ind == 2)   fprintf(file, ['settings.dims = [10];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [2];\n']);     end;
        if (ind == 3)   fprintf(file, ['settings.dims = [10];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [3];\n']);       end;
        if (ind == 4)   fprintf(file, ['settings.dims = [10];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [4];\n']);     end;
        if (ind == 5)   fprintf(file, ['settings.dims = [10];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [5];\n']);     end;
        if (ind == 6)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [3];\n']);  fprintf(file, ['settings.instances = [1];\n']);       end;
        if (ind == 7)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [3];\n']);  fprintf(file, ['settings.instances = [2];\n']);     end;
        if (ind == 8)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [3];\n']);  fprintf(file, ['settings.instances = [3];\n']);     end;
        if (ind == 9)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [3];\n']);  fprintf(file, ['settings.instances = [4];\n']);       end;
        if (ind == 10)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [3];\n']);  fprintf(file, ['settings.instances = [5];\n']);     end;
        if (ind == 11)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [1];\n']);     end;
        if (ind == 12)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [2];\n']);       end;  
      	if (ind == 13)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [3];\n']);     end;
        if (ind == 14)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [4];\n']);     end;
        if (ind == 15)   fprintf(file, ['settings.dims = [20];\n']);   fprintf(file, ['settings.funs = [4];\n']);  fprintf(file, ['settings.instances = [5];\n']);       end;
    end;
    
    if (part == 3)     
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [10,20];\n']);                  
        fprintf(file, ['settings.funs = [3,20,21,22];\n']);      
    end;
    if (part == 4)   
        fprintf(file, ['settings.dims = [10,20];\n']);
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.funs = [23,24];\n']);
    end;
    if (part == 5)
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [2,3,5];\n']);                  
        fprintf(file, ['settings.funs = [1:24];\n']); 
    end;
    if (part == 6)
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [40];\n']);                  
        fprintf(file, ['settings.funs = [1,2,5:14,3];\n']); 
    end;
    if (part == 7)
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [40];\n']);                  
        fprintf(file, ['settings.funs = [15,16,17,18,20];\n']); 
    end;
    if (part == 8)
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [40];\n']);                  
        fprintf(file, ['settings.funs = [19,4];\n']); 
    end;
    if (part == 9)
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [40];\n']);                  
        fprintf(file, ['settings.funs = [21,24];\n']); 
    end;
    if (part == 10)
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [40];\n']);                  
        fprintf(file, ['settings.funs = [22];\n']); 
    end;
    if (part == 11)
        fprintf(file, ['settings.instances = [' inststr '];\n']);
        fprintf(file, ['settings.dims = [40];\n']);                  
        fprintf(file, ['settings.funs = [23];\n']); 
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
    if (part < 6)           fprintf(file, ['settings.MaxEvals = ''1e6*dim'';\n']);
    else                    fprintf(file, ['settings.MaxEvals = ''1e5*40'';\n']);      end; % 1e5 for 40-D
    fprintf(file, ['settings.MaxEvalsWithSurrogate = ''1e5*40'';\n']);
    fprintf(file, ['settings.lambdaMult = ' num2str(lambdaMult) ';\n']);
    fprintf(file, ['settings.muMult = ' num2str(muMult) ';\n']);

    
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