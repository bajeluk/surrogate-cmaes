function [err] = xacmes_estimateModelError(model)

    npt = model.nCrossValidation;
    if (npt == 0)
        err = 0;
        return;
    end;

    xcross = model.CrossValidX;
    fcross = model.CrossValidF;
    
    [arfitness, arindex] = sort(fcross,2,'ascend');
    xcross(:,:) = xcross(:,arindex);
    
    x_test = xcross; 
    npoints = size(x_test,2);
    surrF = evaluateModel(x_test, npoints, model);
       
 	[arfitness, rankarr] = sort(surrF,2,'ascend');
    npt = size(rankarr,2);
    err = EstimError(rankarr,npt);
    if (0) % too slow for npt>>1000
        err = 0;
        for i=1:npt
            cur_err = sum( rankarr(i) > rankarr(i+1:npt) );
            err = err + cur_err;
            if (0) % n^2.0 slow code
                for j=i+1:npt
                    if ( rankarr(i) > rankarr(j))
                        err = err + 1;
                    end;
                end;
            end;
        end;
        err = err / ( npt * (npt - 1) / 2 );
    end;
    

    
