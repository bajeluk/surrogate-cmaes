function Adapter()

global settings;
   
if (settings.CMAactive == 0) && (settings.withSurr == 0) && (settings.withModelOptimization == 0)       name = 'CMA';   end;
if (settings.CMAactive == 0) && (settings.withSurr == 1) && (settings.withModelOptimization == 0)       name = 'ACM';   end;
if (settings.CMAactive == 0) && (settings.withSurr == 0) && (settings.withModelOptimization == 1)       name = 'OptCMA';end;
if (settings.CMAactive == 0) && (settings.withSurr == 1) && (settings.withModelOptimization == 1)       name = 'OptACM';end;
if (settings.CMAactive == 1) && (settings.withSurr == 0) && (settings.withModelOptimization == 0)       name = 'ActiveCMA'; end;
if (settings.CMAactive == 1) && (settings.withSurr == 1) && (settings.withModelOptimization == 0)       name = 'ActiveACM'; end;
if (settings.CMAactive == 1) && (settings.withSurr == 0) && (settings.withModelOptimization == 1)       name = 'ActiveOptCMA';  end;
if (settings.CMAactive == 1) && (settings.withSurr == 1) && (settings.withModelOptimization == 1)       name = 'ActiveOptACM';  end;
datapath = [name settings.pathname];
opt.algName = [name settings.algname];

%experiment
opt.comments = 'COMMENTS';
mkdir(datapath);
gfile = fopen([ './' datapath '/' 'aACM.txt'],'w');
if (gfile == -1)
  disp('gfile is not open');
end;

ggfile = fopen([ './' datapath '/' settings.savfile],'w');
fprintf(ggfile, '\n');
fclose(ggfile);



goalf = 1;
sumspeedup = 0;
nfun = 0;
%tic
name = 'func';

for ifun = settings.funs
    nfun = nfun + 1;
    if (ifun == 1)     	name = 'F1-Sphere';                 end;
    if (ifun == 2)   	name = 'F2-EllipsOrig';             end;
    if (ifun == 3)   	name = 'F3-RastriginOrig';          end;
    if (ifun == 4)   	name = 'F4-RastriginBueche';        end;
    if (ifun == 5)  	name = 'F5-Linear';                 end;
    if (ifun == 6)    	name = 'F6-AttrSector';             end;
    if (ifun == 7)   	name = 'F7-Step-Ellipsoid';         end;
    if (ifun == 8)   	name = 'F8-RosenOrig';              end;
    if (ifun == 9)   	name = 'F9-RosenRot';               end;
    if (ifun == 10)   	name = 'F10-EllipsRot';             end;
    if (ifun == 11) 	name = 'F11-Discus';                end;
    if (ifun == 12)   	name = 'F12-BentCigar';             end;
    if (ifun == 13)  	name = 'F13-SharpRidge';            end;
    if (ifun == 14)    	name = 'F14-SumOfDiffPow';          end;
    if (ifun == 15)   	name = 'F15-RastriginRotated';      end;
    if (ifun == 16)   	name = 'F16-Weierstrass';           end;
    if (ifun == 17)    	name = 'F17-SchafferCond10';        end;
    if (ifun == 18)    	name = 'F18-SchafferCond1000';      end;
    if (ifun == 19)    	name = 'F19-Griewank-Rosenbrock F8F2';	end;
    if (ifun == 20)   	name = 'F20-Schwefelxsinx';         end;
    if (ifun == 21)   	name = 'F21-Gallagher101peaks';     end;
    if (ifun == 22)   	name = 'F22-Gallagher21peaks';      end;
    if (ifun == 23)    	name = 'F23-Katsuuras';             end;
    if (ifun == 24)   	name = 'F24-Lunacek-bi-Rastrigin'; 	end;

ntest2 = 1;

name0 = name;
for ntr=settings.ntarray
%settings.ntr_cur = ntr;
disp(['ntr' num2str(ntr) ]);
%settings.iStep = ntr;
param = zeros(1,1);
for itest1 = settings.dims
 for itest2 = 1:ntest2
    dim = itest1;
    %MAX_EVAL = 1000000*dim;
    MAX_EVAL = myeval(settings.MaxEvals);
    name = [ name0 'D' num2str(dim)];
    avrfevals = 0;
    nsucc = 0;  nunsucc = 0;

    for iinstance = settings.instances
        fgeneric('initialize', ifun, iinstance, datapath, opt); 

        gfile_state_name = [ './' datapath '/' num2str(ifun) '_D' num2str(dim) '_inst' num2str(iinstance) '.txt'];
        settings.gfile_state = fopen(gfile_state_name,'w');
        
        irest = 0;
        settings.iglobalrun = 0;
        while (fgeneric('evaluations') < MAX_EVAL) && (fgeneric('fbest') - fgeneric('ftarget') > 0)
            settings.iglobalrun  = settings.iglobalrun + 1;
            maxeval_available = MAX_EVAL - fgeneric('evaluations');
            %%%%%%  Run x-ACM
            xacmes('MyFunc',dim,maxeval_available);
            %%%%%%
            irest = irest + 1;
            ff = fgeneric('fbest') - fgeneric('ftarget');
        end;
        if (irest > 1)
            disp( [' nrestart:' num2str(irest) ] );
        end;
        fclose(settings.gfile_state);

        if (0)
            disp(sprintf(['  f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                        ' fbest-ftarget=%.4e, elapsed time [h]: %.2f'], ...
                       ifun, dim, iinstance, ...
                       fgeneric('evaluations'), ...
                       fgeneric('fbest') - fgeneric('ftarget')));
        end;
        
        df = fgeneric('fbest') - fgeneric('ftarget');
        if (df < 1e-10)
            nsucc = nsucc + 1;
            avrfevals = avrfevals + fgeneric('evaluations');
        else
            nunsucc = nunsucc + 1;
        end;
        fgeneric('finalize');
    end;
    
    if (nsucc > 0)
        avrfevals = avrfevals / nsucc;
        avrfevals = avrfevals / ( nsucc/(nsucc + nunsucc) );
    else
        avrfevals = 15 * MAX_EVAL;
    end;
    coef = nsucc/(nsucc + nunsucc);
 %   fprintf(gfile, [num2str(param(1)) '\t' num2str(avrfevals) '\n']);
    param(1) = dim;
    fprintf(gfile, [name '\t' num2str(param(1)) '\t' num2str(avrfevals) '\t' num2str(coef) '\n']);
    
    ggfile = fopen([ './' datapath '/' settings.savfile],'a');
    fprintf(ggfile, ['f' num2str(ifun) '\t' num2str(ntr) '\t' num2str(avrfevals) '\n']);
    fclose(ggfile);
    
    speedup = goalf / avrfevals;
    sumspeedup = sumspeedup + speedup;
    disp( [ 'fun:' num2str(ifun) ' dim:' num2str(dim) ' fevals:' num2str(avrfevals) ' nsucc:' num2str(nsucc) ' speedup:'  num2str(speedup) ] );

 end;
% z = toc
end;    
end;    
    if (0)
        SetStyle();
        plot(xx,yy);
    end;
end;

fclose(gfile);

end



