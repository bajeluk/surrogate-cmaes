function client_full_example

    % This program is a full demonstration of a client program
    % for participation in the black box optimization competition.
    % It implements the Nelder Mead algorithm, see
    % http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
    %
    % The algorithm is applied to all problems in a given track;
    % it is demonstrated with the demo account on the trial track.
    % The program shows one possible way to make use of the provided
    % recovery mechanisms. In case of a crash or a power outage the
    % program can be restarted and it will continue the optimization
    % from the saved state. It handles network problems gracefully.

    % configuration
    global LIBNAME; global LIBHFILE;    
    global TRACKNAME;   global USERNAME; global PASSWORD; 
    global LOGIN_DELAY_SECONDS, global LOCK_DELAY_SECONDS;
    LIBNAME = 'bbcomp'; % working directory should include library file(s) and header file
    LIBHFILE = 'bbcomplib.h';   
    LOGFILEPATH = 'logs/';  % path to log files, make sure that it exists and is writable!
    USERNAME = 'demoaccount';
    PASSWORD = 'demopassword';
    TRACKNAME = 'trial';
    LOGIN_DELAY_SECONDS = 10;
    LOCK_DELAY_SECONDS = 60;

    loadBBCompLibrary();

    disp('---------------------------------------------------');
    disp('black box full example competition client in Matlab');
    disp('---------------------------------------------------');
    disp(' ');

    %set configuration options (this is optional)
    result = calllib(LIBNAME,'configure', 1, LOGFILEPATH);
    if (result == 0) 	whenFailed('configure', LIBNAME);   return;	end; 

    safeLogin();

    printTrackNames();   % request the tracks available to this user (this is optional)

    safeSetTrack();

    numProblems = safeGetNumberOfProblems();   % obtain number of problems in the track

    % solve all problems in the track
    for i=0:(numProblems-1)
        solveProblem(i);
    end;

    unloadlibrary(LIBNAME);

function solveProblem(problemID)
    global LIBNAME;
    
    % set the problem
    safeSetProblem(problemID);

    % obtain problem properties
    bud = safeGetBudget(problemID);
    dim = safeGetDimension(problemID);
    evals = safeGetEvaluations(problemID);

    % output status
    if (evals == bud)
        disp(['problem ' num2str(problemID) ' already solved']);
        return;
    else if (evals == 0)
            disp(['problem ' num2str(problemID) ': starting from scratch']);
      	else
            disp(['problem ' num2str(problemID) ': continuing from evaluation ' num2str(evals)]);
      	end;
    end;

    % reserve memory for the simplex and for a number of intermediate points
    simplex = zeros(dim+1,dim+1);
    for j=1:(dim+1)
        for i=1:dim
            % simplex(j,1) = function value
            % simplex(j,2:(dim+1)) = coordinates
            if (i == j) simplex(j,1+i) = 0.8;
            else        simplex(j,1+i) = 0.2;   end;    
        end;
    end;
    simplex_evaluated = 0;
    x0 = zeros(1,dim);
    xr = zeros(1,dim);
    xe = zeros(1,dim);
    xc = zeros(1,dim);

    % recover algorithm state from saved evaluations
    % NOTE: if your algorithm uses random numbers, you may want to take care
    % of seed (see http://fr.mathworks.com/help/matlab/ref/rng.html) 
    % to recover the state of the algorithm exactly
    for e=1:evals
        value = libpointer('doublePtr',1e+100);
        x0 = libpointer('doublePtr',zeros(1,dim));
        [result, x0, value] = calllib(LIBNAME,'history',e-1,x0,value);
        if (result == 0)
            disp('WARNING: history failed.\n'); % note: this evaluation is lost
        else
            if (simplex_evaluated <= dim)
                simplex_evaluated = simplex_evaluated + 1;
                simplex(simplex_evaluated,1) = value;
                for i=1:dim
                     simplex(simplex_evaluated,1+i) = x0(i);
                end;                
            else
                fvals = simplex(1:simplex_evaluated,1);
                [fvals_srt iarr] = sort(fvals,'ascend');
                simplex(1:simplex_evaluated,:) = simplex(iarr,:);
                if (value < simplex(dim+1,1))
                    for i=1:dim
                        simplex(dim+1,1+i) = x0(i);
                    end;
                end;
            end;
        end;
    end;

    % optimization loop
    while (safeGetEvaluations(problemID) < bud)
        if (simplex_evaluated <= dim)
            % evaluate the initial simplex
            simplex_evaluated = simplex_evaluated + 1;
            value = safeEvaluate(problemID, simplex(simplex_evaluated,2:end));
            evals = evals + 1;
            simplex(simplex_evaluated,1) = value;
        else
            % step of the simplex algorithm
            fvals = simplex(1:simplex_evaluated,1);
            [fvals_srt iarr] = sort(fvals,'ascend');
            simplex(1:simplex_evaluated,:) = simplex(iarr,:);
            best = simplex(1,:);
            worst = simplex(end,:);
            disp([num2str(evals) ':' num2str(best(1))]);
            % compute centroid
            for i=1:dim
                x0(1,i) = mean(simplex(1:dim,1+i));   % not 1:(dim+1) ?
            end;

            % reflection
            for i=1:dim
                xr(1,i) = truncate2bounds(2.0 * x0(i) - worst(i+1));
            end;
            f_xr = safeEvaluate(problemID, xr);
            evals = evals + 1;
            if (best(1) <= f_xr && f_xr < worst(1))
                % replace worst point with reflected point
                simplex(end,1) = f_xr;
                simplex(end,2:end) = xr(:);
            else if (f_xr < best(1))
                if (safeGetEvaluations(problemID) >= bud) return;    end;

                % expansion
                for i=1:dim
                    xe(i) = truncate2bounds(3.0 * x0(i) - 2.0 * worst(i+1));
                end;
                f_xe = safeEvaluate(problemID, xe);
                evals = evals + 1;
                if (f_xe < f_xr)
                    % replace worst point with expanded point
                    simplex(end,1) = f_xe;
                    simplex(end,2:end) = xe(:);
                else
                    % replace worst point with reflected point
                    simplex(end,1) = f_xr;
                    simplex(end,2:end) = xr(:);
                end;
            else   
                if (safeGetEvaluations(problemID) >= bud) return;    end;

                % contraction
                for i=1:dim
                    xc(i) = truncate2bounds(0.5 * x0(i) + 0.5 * worst(i+1));
                end;
                f_xc = safeEvaluate(problemID, xc);
                evals = evals + 1;
                if (f_xc < worst(1))
                    % replace worst point with contracted point
                    simplex(end,1) = f_xc;
                    simplex(end,2:end) = xc(:);
                else
                    % reduction
                    for j=1:(dim+1)
                        if (safeGetEvaluations(problemID) >= bud) return;   end;
                        for i=1:dim
                            simplex(j,i+1)= truncate2bounds(0.5 * best(i+1) + 0.5 * simplex(j,i+1));
                        end;
                        simplex(j,1) = safeEvaluate(problemID, simplex(j,2:end));
                        evals = evals + 1;
                    end;
                end;
            end;
            end;
        end;
    end;
    
    
%%% some auxiliary functions %%%    
    
function whenFailed(funcname, LIBNAME)
    disp([funcname '() failed: ' calllib(LIBNAME,'errorMessage')]);    
    unloadlibrary(LIBNAME);

function loadBBCompLibrary()
    global LIBNAME; global LIBHFILE;
    if (libisloaded(LIBNAME) == 1) % first unload if already loaded
        unloadlibrary(LIBNAME);
    end;
    loadlibrary(LIBNAME,LIBHFILE);
    Library_Functions = libfunctions(LIBNAME, '-full') % print library functions
    if (numel(Library_Functions) == 0)
        disp('error loading library');
        return;
    end;

    
function safeLogin()
    global LIBNAME; global USERNAME; global PASSWORD; global LOGIN_DELAY_SECONDS;
    while (1)
        result = calllib(LIBNAME,'login',USERNAME,PASSWORD); 
        if (result ~= 0) return;    end;
        msg = calllib(LIBNAME,'errorMessage');
        disp(['WARNING: login failed: ' msg]);
        pause(LOGIN_DELAY_SECONDS);
        if (strcmp(msg, 'already logged in') == 0) return;  end;
    end

    
function safeSetTrack()
    global TRACKNAME; global LIBNAME; global USERNAME; 
    global PASSWORD; global LOGIN_DELAY_SECONDS;
    while (1)
        result = calllib(LIBNAME,'setTrack',TRACKNAME);
        if (result ~= 0) return;    end;
        disp(['WARNING: setTrack failed: ' calllib(LIBNAME,'errorMessage')]);
        safeLogin();
    end;

    
function printTrackNames()
    global LIBNAME;
    numTracks = calllib(LIBNAME,'numberOfTracks');     
    if (numTracks == 0) whenFailed('numberOfTracks', LIBNAME);    return;	end;
    disp([num2str(numTracks) ' track(s):']);
    for i=1:numTracks
        trackname = calllib(LIBNAME,'trackName',i-1); % indexing from 0
        if (isempty(trackname))         whenFailed('trackName', LIBNAME);       return;
        else                            disp([' ' num2str(i) ': ' trackname]);	end;
    end;

    
function numProblems = safeGetNumberOfProblems()
    global LIBNAME;
    numProblems = calllib(LIBNAME,'numberOfProblems');
    if (numProblems == 0) whenFailed('numberOfProblems', LIBNAME);     return;	end;
    disp(['The track consists of ' num2str(numProblems) ' problems.']);
  
    
function safeSetProblem(problemID)
    global LIBNAME; global LOCK_DELAY_SECONDS;
    while (1)
        result = calllib(LIBNAME,'setProblem',problemID);
        if (result ~= 0) return;    end;
        msg = calllib(LIBNAME,'errorMessage');
        disp(['WARNING: setProblem failed: ' msg]);
        %if (strcmp(msg, 'failed to acquire lock') == 0) 
            pause(LOCK_DELAY_SECONDS);
       % else
            safeSetTrack();
       % end; 
    end;
  
    
function result = safeGetDimension(problemID)
    global LIBNAME;
    while (1)
        result = calllib(LIBNAME,'dimension');
        if (result ~= 0) return;    end;
        disp(['WARNING: dimension failed: ' calllib(LIBNAME,'errorMessage')]);
        safeSetProblem(problemID);
    end;
   
    
function result = safeGetBudget(problemID)
    global LIBNAME;
    while (1)
        result = calllib(LIBNAME,'budget');
        if (result ~= 0) return; end;
        disp(['WARNING: budget failed: ' calllib(LIBNAME,'errorMessage')]);
        safeSetProblem(problemID);
    end;
   
    
function result = safeGetEvaluations(problemID)
    global LIBNAME;
    while (1)
        result = calllib(LIBNAME,'evaluations');
        if (result >= 0) return;    end;
        disp(['WARNING: evaluations failed: ' calllib(LIBNAME,'errorMessage')]);
        safeSetProblem(problemID);
    end;
    
function value = safeEvaluate(problemID, point)
    global LIBNAME;
	while (1)
        % query the black box
        value = libpointer('doublePtr',1e+100);
        [result, TMP, value] = calllib(LIBNAME,'evaluate',point,value);
		if (result ~= 0) return;    end;
		disp(['WARNING: evaluate failed: ' calllib(LIBNAME,'errorMessage')]);
		safeSetProblem(problemID);
    end;
    
function value = truncate2bounds(value)
    value = max(0.0, min(1.0, value));
    