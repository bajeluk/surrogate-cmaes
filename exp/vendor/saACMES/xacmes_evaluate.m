function Fit = xacmes_evaluate(x)

global algo;
global settings;

npoints = size(x,2);
if (algo.realfunc == 1)     % using Real Fitness Function F(x)
    Fit = feval(settings.strfitnessfct, x);
    
    
    if (settings.withSurr == 1)
        for i=1:npoints
            algo.nevals = algo.nevals + 1;
            algo.iSZ = algo.iSZ + 1;
            algo.aSZ = algo.aSZ + 1;
            if (algo.aSZ > algo.maxArchSize)    
                algo.aSZ = algo.maxArchSize;
            end;
            if (algo.iSZ > algo.maxArchSize)    
                algo.iSZ = 1;
                if (0)
                    algo.iSZ = algo.aSZ;
                    algo.ARX(:,1) = [];
                    algo.ARF(1) = [];
                end;
            end

            algo.ARX(:,algo.iSZ) = x(:,i);  
            algo.ARF(algo.iSZ) = Fit(i);
            algo.ARneval(algo.iSZ) = algo.nevals;
        end;
    else
        algo.nevals = algo.nevals + npoints;
    end;
    [minFit,iminX] = min(Fit);
    if (minFit < algo.Fmin)
        algo.Fmin = minFit;
        algo.xmin = x(:,iminX);
    end;
end;   
if (algo.realfunc == 0)     % using Surrogate Fitness Function F^(x)
    global curModel;

 	x_test = x;
%    Fit = evaluateModel(x_test, npoints, curModel);
    global Models;
    if ((Models.imodel == 1) || (npoints == 1))
        Fit = evaluateModel(x_test, npoints, curModel);
    end;
    if ((Models.imodel > 1) && (npoints > 1))
        Fit = zeros(1,npoints);
        for i=1:Models.imodel
            locFit = evaluateModel(x_test, npoints, Models.models(i));
            [locFitsort, idx] = sort(locFit);
            locFitrank = 1:npoints;
            locFitrank(idx) = 1:npoints;
            Fit = Fit + Models.weights(i) * locFitrank;
        end;
        zzzzz = 0;
        %%Fit = evaluateModel(x_test, npoints, curModel);
    end;
end;
if (algo.realfunc == 2)     % using Model Quality Function H(theta)
    global optModel;
    global algo;
    Fit = 0;
    nx = size(x,1);
    penalty = 0;
    correct = 1;
    for i=1:nx
        if (x(i) < 0)
            penalty = penalty + abs(0 - x(i));
        end;
        if (x(i) > 1)
            penalty = penalty + abs(1 - x(i));
        end;
        if (x(i) < 0.0) || (x(i) > 1.0)
            correct = 0;
        end;
    end;
    
    for i=1:nx
        maxV = optModel.xmax(i);
        minV = optModel.xmin(i);
        if (i == 1)
            if (maxV > algo.aSZ)
      %          max = algo.aSZ;
            end;
        end;
        coeff(i) = minV + x(i,1) * (maxV - minV);
    end;
    if (coeff(1) > algo.aSZ);
        penalty = coeff(1) - (algo.aSZ);
        correct = 0;
    end;
    
    if (correct == 0)
        Fit = 1 + penalty;
        return;
    end;

    ModelType = 1;
    model = buildModel(optModel.algo_initial_state, coeff, ModelType);
    
    %estimate f^(x) values 
    if (model.nCrossValidation == 0)
        model.CrossValidX = optModel.CrossValidX;
        model.CrossValidF = optModel.CrossValidF;
    else
        for i=1:optModel.nCrossValidation
            model.CrossValidX(:,model.nCrossValidation + i) = optModel.CrossValidX(:,i);
            model.CrossValidF(model.nCrossValidation + i) = optModel.CrossValidF(i);
     	end;
    end;
    model.nCrossValidation = model.nCrossValidation + optModel.nCrossValidation;
    
    err = xacmes_estimateModelError(model);
    
    Fit = err;
    
end;

