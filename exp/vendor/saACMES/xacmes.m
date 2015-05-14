function [xbest, y_eval] = xacmes(strfitnessfct, N, MAX_EVAL)

global settings;
settings.strfitnessfct = strfitnessfct;
settings.ftarget = fgeneric('ftarget');
wrapperfct = 'xacmes_evaluate';

opts = cmaes_initialize('defaults');
opts.MaxFunEvals = MAX_EVAL;
opts.StopFitness = 0;
opts.DispModulo = -1;
opts.DispFinal = 'off';
opts.EvalParallel       = 'yes  % original: no   % objective function FUN accepts NxM matrix, with M>1?';

if (settings.BIPOP == 1) % as suggested for BIPOP on BBOB-2009
   if (settings.noisy == 0)
        opts.MaxIter   	= '100 + 50 * (N+3)^2 / sqrt(popsize)  % original: 1e3*(N+5)^2/sqrt(popsize) % maximal number of iterations per (re)start';
        opts.StopOnStagnation	= 'on  % stop when fitness stagnates for a long time';
   end;
   opts.TolFun             = '1e-12 % stop if fun-changes smaller TolFun';
   opts.TolHistFun         = '1e-12 % stop if back fun-changes smaller TolHistFun';
   opts.TolX               = '1e-12*max(insigma) % stop if x-change smaller TolX'; % TODO: re-consider, for sharpR 1e-8 is too tight, dito 1e-9 for ackley and ftarget=1e-8, rastrigin with fnoise3 needs > 1e-10
   opts.TolUpX             = 'Inf  % original: 1e3*max(insigma) % stop if x-changes larger TolUpX';
   opts.CMA.ccov1   = 'min(2,lambda/3) / ((N+1.3)^2+mueff)  % learning rate for rank-one update'; 
   opts.CMA.ccovmu  = 'min(2,lambda/3) * (mueff-2+1/mueff) / ((N+2)^2+mueff) % learning rate for rank-mu update'; 
   opts.CMA.ccum = '((N+4 + 2*mueff/N) / (4 + mueff/N))^-1 '; 
end;

if (N < 10)
    settings.lambdaMult = 1;
end;
if (N >= 10)
    settings.lambdaMult = floor(10^(1 + log2(N/10)));
end;

opts.CMA.active = settings.CMAactive;
opts.Restarts = 9;
%opts.LogPlot = 'on';
opts.LogModulo = 0;
opts.LogTime   = 0;
opts.SaveVariables = 'off';
opts.ReadSignals = 'off';

sigma0 = 2.0;
if (settings.noisy == 1)
    sigma0 = 2.0;
    noisy_coef = 1/5;
    opts.CMA.ccov1 = [num2str(noisy_coef) '*' opts.CMA.ccov1];
    opts.CMA.ccovmu = [num2str(noisy_coef) '*' opts.CMA.ccovmu];
    opts.TolX = '0*max(insigma) % stop if x-change smaller TolX';
end;
%opts.PopSize = 6;
    
xstart = [ '-4 + 8*rand(' num2str(N) ',1)' ];    
cur_state = cmaes_initialize(wrapperfct, xstart, sigma0, opts);


if (settings.withFileDisp == 1)
    if (settings.withModelOptimization == 0) && (settings.withSurr == 1)
        fprintf(settings.gfile_state, [ 'nevals' '\t' 'iter' '\t' 'est_err' '\t' 'Fmin' '\t' 'iSTEP' '\n']);
    end;
    if (settings.withModelOptimization == 1) && (settings.withSurr == 1)
        fprintf(settings.gfile_state, [ 'nevals' '\t' 'iter' '\t' 'est_err' '\t' 'Fmin' '\t' 'iSTEP' '\t' ...
            'coeffs(:,1)' '\t' 'coeffs(:,2)' '\t' 'coeffs(:,3)' '\n']);
    end;               
end;


%%%**** for BIPOP ****%%%%
if (settings.BIPOP == 1)
    nrunswithsmallpopsize = 1;
    budget.smallpopsi = [];
    budget.largepopsi = [];
end;
SmallPOPdivBigPOPbudget = 1.0;
%%%**** for BIPOP ****%%%%

nrestarts = myeval(cur_state.opts.Restarts);
bipop_criterion = 0;

%%%**** for new restart BIPOP/IPOP ****%%%%
if (settings.newRestartRules == 1)
    SmallPOPdivBigPOPbudget = 1.0;
    n_resultsBigPOP = 0;
    n_resultsSmallPOP = 0;
    resultsBigPOP = zeros(1,1);
    resultsSmallPOP = zeros(1,1);
end;
%%%**** for new restart BIPOP/IPOP ****%%%%


y_eval = [];  % BAJELUK BEST/COUNTEVAL RECORDING


while cur_state.irun <= nrestarts || bipop_criterion % loop with restarts 
    %%% set default parameters     .begin
    global zero1;   global model;   model = zero1;
    global zero2;   global algo;    algo = zero2;
    global zero3;   global curModel;    curModel = zero3;
    global zero4;   global optModel;  optModel = zero4;
    global zero5;   global Models;  Models = zero5;

    algo.nevals = 0;
    algo.realfunc = 1;
    algo.iSZ = 0;
    algo.aSZ = 0;
    algo.Fmin = 1e+30;
     
    modelType = settings.modelType;
    optModel = InitModelParameters(modelType,N,settings.withModelOptimization);
        
    algo.CMAactive = settings.CMAactive;
    algo.withDisp = settings.withDisp;
    algo.withSurr = settings.withSurr;
    algo.withModelOptimization = settings.withModelOptimization;
    algo.iterstart = settings.iterstart;
    algo.withFileDisp = settings.withFileDisp;

    algo.maxArchSize = optModel.MaxTrainingPoints;
    algo.ARX = zeros(N, algo.maxArchSize);  
    algo.ARF = zeros(1, algo.maxArchSize); 
    algo.ARneval = zeros(1, algo.maxArchSize);
    
    
    cur_state.irun = cur_state.irun + 1; 
    if (settings.newRestartRules == 1) 
        if ((cur_state.irun > 3) || (settings.iglobalrun > 1)) %% do not use surrogate for multirestarts ( inrun > 2)
            algo.withSurr = 0;
            algo.withModelOptimization = 0;
        end;
    end;

    
    cur_state = cmaes_initializeRun(cur_state);
    if (settings.BIPOP == 1)
        lambda0 = floor(myeval(opts.PopSize) * myeval(opts.IncPopSize)^(cur_state.irun-nrunswithsmallpopsize));
        popsize = lambda0;
        cur_state.lambda = lambda0;
    	cur_state.popsize = popsize;
    end;
    cur_state.stopflag = {};
    
    algo.err = 0.5;
    algo.iter = 0; 
    iSTEP = 0;
    WeDontBelieveInModel_STEP = 0;
    
    isIPOPRun = 1;
   	%%%**** for BIPOP ****%%%%
    if (settings.BIPOP == 1)
          insigmafac = 1;
          % restarts with small popsize, needs to be refined and cleaned
        %  if 1 < 3 && cur_state.irun > 2 && 1 * sum(budget.smallpopsi) < sum(budget.largepopsi) 
          if 1 < 3 && cur_state.irun > 2 && sum(budget.smallpopsi) < SmallPOPdivBigPOPbudget*sum(budget.largepopsi)     %%%**** for new restart BIPOP/IPOP ****%%%%
            nrestarts = nrestarts + 1; 
            nrunswithsmallpopsize = nrunswithsmallpopsize + 1; 
            budget.irunwithsmallpopsize = cur_state.irun; 
                   %%%%%% qqqqqqqTODO  \/
            insigmafac = 0.01^(rand(1)^1);  % insigma should not be changed

            % set lambda0
            lambda0 = myeval(opts.PopSize);
            if (settings.newRestartRules == 0)
                lambda0 = floor(lambda0 * (cur_state.lambda/myeval(opts.IncPopSize)/lambda0)^(rand(1)^2));
            end;
            %%%**** for new restart BIPOP/IPOP ****%%%%
            if (settings.newRestartRules == 1)
                lambda0 = myeval(opts.PopSize);  % restart with default population size, but randomized sigma 
            end;
            %%%**** for new restart BIPOP/IPOP ****%%%%
            popsize = lambda0;
            lambda = lambda0;
            budget.maxiter = 0.5 * sum(budget.largepopsi) / lambda0; % CAVE: not very precise with adaptive popsize
            cur_state.sigma = insigmafac*max(cur_state.insigma);  % overall standard deviation 
            cur_state.lambda = lambda;
            cur_state.popsize = popsize;
            isIPOPRun = 0;
          else
            budget.maxiter = inf;  % for setting of maxiter
          end
           
          cur_state.stopMaxIter = min(myeval(opts.MaxIter), budget.maxiter);
    end;
    %%%**** for BIPOP ****%%%%
    
    %%%**** for new restart BIPOP/IPOP ****%%%%
    if ((settings.newRestartRules == 1) && (isIPOPRun == 1))
    	if (settings.BIPOP == 0) nrunswithsmallpopsize = 1;     end;
    	ibigrun = cur_state.irun-nrunswithsmallpopsize + 1;
     	insigmafac = 1.6^(-(ibigrun-1));%1.6^(-(ibigrun-1));
     	cur_state.sigma = insigmafac*max(cur_state.insigma);  
    end;
    %%%**** for new restart BIPOP/IPOP ****%%%%
   
    if (cur_state.irun > 1)
        disp( ['restart #' num2str(cur_state.irun-1) ' popsize:' num2str(cur_state.popsize) ' sigma:' num2str(cur_state.sigma)  ] );
    end;
    
    %%% set default parameters     .end

    stop = 0;
    while (stop == 0) %%(algo.nevals < 10000) && (algo.Fmin > 1e-10) && 

        if (fgeneric('evaluations') > myeval(settings.MaxEvalsWithSurrogate))
            algo.withSurr = 0;
        end;
        
        if (algo.withSurr == 0) % original CMA
            cur_state = cmaes_iteration(cur_state);
            stop = cur_state.stop;
        end;
        if (algo.withSurr == 1) % evaluate lambda points on some iSTEP'th step with adaptation of iSTEP
            
            err = 0;
            algo.iter = algo.iter + 1;
            %if (algo.iter > 30)
           %     settings.lambdaMult = floor(algo.iter*0.03);
           %     if (settings.lambdaMult < 1)    settings.lambdaMult = 1;    end;
           %     if (settings.lambdaMult > 10)    settings.lambdaMult = 10;    end;
           % end;
            
            if (algo.iter < algo.iterstart)
                cur_state = cmaes_iteration(cur_state);
                algo.sav_fitness = cur_state.fitness;
                algo.sav_out = cur_state.out;
                algo.sav_counteval = cur_state.counteval;
            else   
                algo.realfunc = 0;

                %%% recover current coefficients/hyper-parameters    .begin
                for i=1:size( optModel.coeff, 2);
                    cur_coeff(i) = (optModel.coeff(i) - optModel.xmin(i)) / ( optModel.xmax(i) - optModel.xmin(i));
                end;
                algo.coeffs(algo.iter,:) = cur_coeff(1,:);
                algo.coeffs_mean = mean( algo.coeffs(algo.iterstart:algo.iter,:) );
                %%% recover current coefficients/hyper-parameters    .end

                %%% BUILD A SURROGATE MODEL     .begin
                modelType = settings.modelType;
                model = buildModel(cur_state,optModel.coeff,modelType);
                global curModel;
                curModel = model;
                    %%%
                Models.imodel = 1;
                if (settings.withModelEnsembles == 1)
                    Models.models(Models.imodel) = model;
                    Models.weights(Models.imodel) = 1;
                        if (optModel.init == 1)
                        if (optModel.cur_state.countiter > 1)

                            for i=1:optModel.cur_state.mu
                                Models.imodel = i;
                                curparams = optModel.cur_state.Xnew_sorted(:,i);                
                                %xmean = optModel.cur_state.Xnew_sorted(:,1) ;    
                                ncoeff = size( optModel.coeff, 2);
                                correct = 1;
                                for j=1:ncoeff
                                    if (curparams(j) < 0) || (curparams(j) > 1)
                                        correct = 0;
                                    end;
                                    pcoeff(j) = optModel.xmin(j) + curparams(j) * ( optModel.xmax(j) - optModel.xmin(j));
                                end;
                                if (correct == 1)
                                    Models.coeffs(i,:) = pcoeff;
                                else
                                    Models.coeffs(i,:) = optModel.def_coeff(:);
                                end; 
                                icurmodel = buildModel(cur_state,Models.coeffs(i,:),modelType);
                                Models.weights(i) = optModel.cur_state.weights(i);
                                Models.models(i) = icurmodel;
                            end;
                        end;
                        end;
                end;
                %%%
                %%% BUILD A SURROGATE MODEL     .end

                initial_state = cur_state;

                %%% OPTIMIZE SURROGATE F^(X) for iSTEP iterations    .begin
                cur_state.savemodulo = 0;
                
                cur_lambdaMult = settings.lambdaMult;
                cur_muMult = settings.muMult;
                if (iSTEP <= settings.largeLambdaMinIter)
                    cur_lambdaMult = 1;
                end;
                
                sav_mu = cur_state.opts.ParentNumber;
                popsize = cur_state.lambda;
                cur_mu = myeval(cur_state.opts.ParentNumber);
                maxpop = max(10*myeval(opts.PopSize) * cur_lambdaMult, 100000);
                cur_state.lambda = cur_state.lambda * cur_lambdaMult;
                cur_state.popsize = cur_state.lambda;
                if (cur_state.popsize > maxpop)  cur_state.popsize = maxpop;  end;
                new_mu = cur_mu * cur_muMult;
                cur_state.opts.ParentNumber = num2str(new_mu);
              %  iSTEP = settings.iStep;
                for i=1:iSTEP
                    algo.realfunc = 0;
                    cur_state = cmaes_iteration(cur_state);
                end;
                cur_state.lambda = cur_state.lambda / cur_lambdaMult;
                cur_state.popsize = cur_state.lambda;
                cur_state.opts.ParentNumber = sav_mu;
                %%% OPTIMIZE F(X) for iSTEP-1 iterations    .end
                

                algo_without_new_points = algo;

                %%% ONE GENERATION ON F(X)      .begin
                algo.realfunc = 1;
                cur_state.savemodulo = 1;
                cur_state.fitness = algo.sav_fitness;
                cur_state.out = algo.sav_out;
                cur_state.counteval = algo.sav_counteval;

                cur_state = cmaes_iteration(cur_state);

                algo.sav_fitness = cur_state.fitness;
                algo.sav_out = cur_state.out;
                algo.sav_counteval = cur_state.counteval;
                %%% ONE GENERATION ON F(X)      .end

                stop = cur_state.stop;
                algo_with_new_points = algo;

                %%% Estimate Model Error        .begin
                xnew = cur_state.Xnew_sorted;
                fnew = cur_state.arfitness;
                npt = size( xnew, 2);

                if (model.nCrossValidation == 0)
                    model.CrossValidX = xnew;
                    model.CrossValidF = fnew;
                else
                    for i=1:npt
                        model.CrossValidX(:,model.nCrossValidation + i) = xnew(:,i);
                        model.CrossValidF(model.nCrossValidation + i) = fnew(i);
                    end;
                end;
                model.nCrossValidation = model.nCrossValidation + npt;

                if (Models.imodel == 1)
                    err = xacmes_estimateModelError(model);
                else
                    sav = algo.realfunc;
                    algo.realfunc = 0;
                    
                     err = xacmes_estimateModelError(model);
                    
                    algo.realfunc = sav;
                end;
                algo.model_err(algo.iter) = err;
                %%% Estimate Model Error        .end
                
                %%% If iSTEP < 1 we suppose that the model is "bad" and  
                %%% we don't believe hyper-parameter optimization can improve it, 
                %%% so we don't use it for hyper_lambda iterations.
                if (((iSTEP < settings.iSTEPminForHyperOptimization) && (WeDontBelieveInModel_STEP < settings.hyper_lambda) && (cur_state.irun > 1)) ) 
                    WeDontBelieveInModel = 1;
                    WeDontBelieveInModel_STEP = WeDontBelieveInModel_STEP + 1;
                else
                    WeDontBelieveInModel = 0;
                    WeDontBelieveInModel_STEP = 0;
                end;
                
                
                if ((algo.withModelOptimization == 1)  && (WeDontBelieveInModel == 0))
                    algo = algo_without_new_points; %otherwise the training set will include the test and models will be almost ideal

                %	optModel.x = cur_state.Xnew_sorted(:,:);
                    optModel.nCrossValidation = npt;
                    optModel.CrossValidX = xnew;
                    optModel.CrossValidF = fnew;
                    optModel.nx = cur_state.lambda;
                    optModel.algo_initial_state = initial_state;
                    algo.realfunc = 2;
                %  	if (rem(algo.iter,40) == 0) % to restart hyper-parameter optimization 
                %        optModel.init = 0;
                %  	end;
                    if (optModel.init == 0)   %Intialization of CMA-ES for hyper-parameter optimization
                        hyper_N = size( optModel.coeff, 2);
                        hyper_X_a = 0.2;   hyper_X_b = 0.8;    hyper_sigma = 0.6;    hyper_lambda = settings.hyper_lambda;   

                        hyper_opts = cmaes_initialize('defaults');
                        hyper_opts.MaxFunEvals = 1e+10;
                        hyper_opts.StopFitness = -1;
                        hyper_opts.DispModulo = -1;
                        hyper_opts.DispFinal = 'off';
                        hyper_opts.CMA.active = algo.CMAactive;
                        %hyper_opts.TolFun       = 1e-15;
                        %hyper_opts.TolHistFun   = 1e-15;
                        hyper_opts.LogModulo = 0;
                        hyper_opts.LogTime   = 0;
                        hyper_opts.SaveVariables = 'off';

                        hyper_opts.PopSize = hyper_lambda;

                        hyper_xstart = [ num2str(hyper_X_a) ' + ' num2str(hyper_X_b-hyper_X_a) '*rand(' num2str(hyper_N) ',1)' ];

                        optModel.cur_state = cmaes_initialize(wrapperfct, hyper_xstart, hyper_sigma, hyper_opts);
                        optModel.cur_state.irun = optModel.cur_state.irun + 1;
                        optModel.cur_state = cmaes_initializeRun(optModel.cur_state); 
                        optModel.cur_state.stopflag = {};
                        optModel.init = 1;
                        optModel.coeff = optModel.def_coeff;
                    else        % One Iteration
                        optModel.cur_state = cmaes_iteration(optModel.cur_state);
                        
                        hyper_stop = cur_state.stop;
                        if (hyper_stop == 1)
                       %     disp('hyper_stop');
                       %     zz = 0;
                        end;

                        xmean = optModel.cur_state.xmean;                 % new model parameters will be chosen as the mean of current distribution
                        %xmean = optModel.cur_state.Xnew_sorted(:,1) ;    % or as the best solution among lambda recently evaluated models
                        ncoeff = size( optModel.coeff, 2);
                        correct = 1;
                        for j=1:ncoeff
                            if (xmean(j) < 0) || (xmean(j) > 1)
                                correct = 0;
                            end;
                            coeff(j) = optModel.xmin(j) + xmean(j) * ( optModel.xmax(j) - optModel.xmin(j));
                        end;
                        if (correct == 1)
                            optModel.coeff = coeff;
                        else
                            optModel.coeff = optModel.def_coeff;
                        end;

                        if (algo.withDisp == 1)
                            str = '';
                            str_mean = '';
                            for j=1:ncoeff
                                str = [ str ' ' num2str( xmean(j,1) ) ];
                                str_mean = [ str_mean ' ' num2str( algo.coeffs_mean(j) ) ];
                            end;
                            arrf = optModel.cur_state.arfitness;
                            disp( [ 'xmean: ' str ] );
                            disp( [ 'coeffs_mean: ' str_mean] );
                            disp( [ 'fit:' num2str(cur_state.arfitness) ] );
                            disp( [ 'sigma:' num2str(optModel.cur_state.sigma)] );
                            disp( [ 'cross:' num2str(arrf) ] );
                        end;
   
                    end;
                    algo = algo_with_new_points;
                    algo.realfunc = 1;
                end;
            end;

            
            if (algo.iter < algo.iterstart)
                algo.model_avrerr(algo.iter) = algo.err;
            end;
            if (algo.iter >= algo.iterstart)

                %%% adapt iSTEP     .begin
                est_err = err;
                algo.err = algo.err*(1-settings.alpha) + settings.alpha*est_err;
                algo.model_avrerr(algo.iter) = algo.err;
              	pos = (settings.maxerr - algo.err) / settings.maxerr;   % using this simple rule
                iSTEP = floor( pos * settings.maxStepts ) - 1;  
                if (iSTEP < 0)                      iSTEP = 0;                      end;
                if (iSTEP > settings.maxStepts)     iSTEP = settings.maxStepts;     end;
                %%% adapt iSTEP     .end
                
                %%% print to Figure/File    .begin
                szcoeffs = size( algo.coeffs,2 );
                if ((algo.withDisp == 1) || (settings.withFileDisp == 1))
                    if (algo.withDisp == 1)
                        disp( [ 'nevals: ' num2str(algo.nevals) ' fit: ' num2str(algo.Fmin) ' model_error: ' num2str(est_err) ' istep: ' num2str(iSTEP) ] );

                        % SetStyle();
                        subplot(2,1,1);
                        curlambda = cur_state.lambda;
                        plot( (1:algo.iter)*curlambda, algo.model_avrerr(:), 'color', 'blue' ); hold off;   
                        xlabel('Number of evaluations');    ylabel('Model error');
                        subplot(2,1,2);
                        colors = {'blue','red','black','green','magenta','cyan'};
                        for i=1:szcoeffs
                            color = colors{i};
                            plot( (1:algo.iter)*curlambda, algo.coeffs(:,i), 'color', color ); hold on;
                        end;
                        xlabel('Number of evaluations');    ylabel('Coefficient value');
                        if (modelType == 1)
                            if (szcoeffs == 3)                legend('nTraining','C0','CPower','Orientation','horizontal'); end;
                            if (szcoeffs == 4)                legend('nTraining','C0','CPower','sigmaA','Orientation','horizontal'); end;
                            if (szcoeffs == 5)                legend('nTraining','C0','CPower','sigmaA','IterCoef','Orientation','horizontal'); end;
                        end;
                        if (modelType == 2)
                            if (szcoeffs == 3)                legend('nTraining','C0','CPower','Orientation','horizontal'); end;
                            if (szcoeffs == 4)                legend('nTraining','C0','CPower','sigmaA','Orientation','horizontal'); end;
                            if (szcoeffs == 5)                legend('nTraining','C0','CPower','sigmaA','PowerCoeff','Orientation','horizontal'); end;
                            if (szcoeffs == 6)                legend('nTraining','C0','CPower','sigmaA','PowerTrans','PowerEps','Orientation','horizontal'); end;
                        end;
                        hold off;
                    end;
                    if (settings.withFileDisp == 1)
                        if (settings.withModelOptimization == 0) && (algo.withSurr == 1)
                            fprintf(settings.gfile_state, [ num2str(algo.nevals) '\t' num2str(algo.iter) '\t' num2str(algo.model_avrerr(algo.iter)) ...
                                '\t' num2str(algo.Fmin) ' ' num2str(iSTEP) '\n']);
                        end;
                        if (settings.withModelOptimization == 1) && (algo.withSurr == 1)
                            str = [ num2str(algo.nevals) '\t' num2str(algo.iter) '\t' num2str(algo.model_avrerr(algo.iter)) '\t' num2str(algo.Fmin) '\t' num2str(iSTEP) '\t' ];
                           for i=1:szcoeffs-1
                               str = [ str num2str(algo.coeffs(algo.iter,i)) '\t' ];
                           end;
                           str = [ str num2str(algo.coeffs(algo.iter,szcoeffs)) '\n' ];
                           fprintf(settings.gfile_state, str);
                        end;                
                    end;
                end;
                %%% print to Figure/File    .end
            end;
        end;

        if (algo.withDisp == 1) && (rand() < 0.1)
            disp( [ num2str(algo.nevals) ' ' num2str( algo.Fmin ) ] );
        end;

        y_eval = [y_eval; algo.Fmin algo.nevals];  % BAJELUK BEST/COUNTEVAL/SURROGATE_STATS RECORDING
    end;

    [cur_state, break_result] = cmaes_finalize(cur_state);
    if (break_result == 1)
        break
    end
    
    
    %%%**** for new restart BIPOP/IPOP ****%%%%
    if (settings.newRestartRules == 1)
        if (isIPOPRun == 1) || (cur_state.irun < 3)
            n_resultsBigPOP = n_resultsBigPOP + 1;
            resultsBigPOP(n_resultsBigPOP) = algo.Fmin;
        end;
        if (isIPOPRun == 0) || (cur_state.irun < 3)
            n_resultsSmallPOP = n_resultsSmallPOP + 1;
            resultsSmallPOP(n_resultsSmallPOP) = algo.Fmin;
        end;
        if (min(resultsBigPOP) < min(resultsSmallPOP))
            SmallPOPdivBigPOPbudget = 0.5;
        else
            SmallPOPdivBigPOPbudget = 2.0;
        end;
    end;
    %%%**** for new restart BIPOP/IPOP ****%%%%
    
    %%%**** for BIPOP ****%%%%
    if (settings.BIPOP == 1)
        % compute budgets for restarts depending on popsize
        if cur_state.irun == 1
          budget.smallpopsi = cur_state.counteval;
          budget.largepopsi = [];
          budget.counteval0 = 0; 
          budget.irunwithsmallpopsize = 1;
        end
        if budget.irunwithsmallpopsize == cur_state.irun
          budget.smallpopsi(end+1) = cur_state.counteval - budget.counteval0; 
        elseif cur_state.irun > 1  % first run does not count in the large budget
          budget.largepopsi(end+1) = cur_state.counteval - budget.counteval0; 
        end
        budget.counteval0 = cur_state.counteval;
        %        bipop_criterion = sum(budget.smallpopsi) < sum(budget.largepopsi);
        bipop_criterion = sum(budget.smallpopsi) < SmallPOPdivBigPOPbudget*sum(budget.largepopsi);  %%%**** for new restart BIPOP/IPOP ****%%%%
    end;
    %%%**** for BIPOP ****%%%%
    if (stop == 1)
  %      disp(['stopflag: ' cur_state.stopflag]);      
    end;
end;

xbest = cur_state.xmean;
