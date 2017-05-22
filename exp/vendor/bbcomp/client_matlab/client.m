function client
	% change the track name to trialMO for multi-objective optimization
	track = 'trial'
	% track = 'trialMO'

    % working directory should include library file(s) and header file
	% for Windows platforms: bbcomp.dll (32 or 64 bit), msvcr120.dll, msvcp120.dll
    libname = 'bbcomp';
    if (libisloaded(libname) == 1)
        unloadlibrary(libname);
    end;
    loadlibrary(libname,'bbcomplib.h', 'mfilename', 'mHeader');
    Library_Functions = libfunctions(libname, '-full') % print library functions

    disp('----------------------------------------------');
    disp('black box example competition client in Matlab');
    disp('----------------------------------------------');
    disp(' ');

    %set configuration options (this is optional)
    result = calllib(libname,'configure', ...
                    1, ...              %   enable history (the default)
                    './logs/' ...   	%   path to log files, make sure that it exists and is writable!
                );
    if (result == 0) 	whenFailed('configure', libname);   return;	end;          

    % login with demo account - this should grant access to the "trial" track (for testing and debugging)
    result = calllib(libname,'login','demoaccount','demopassword');    if (result == 0)    whenFailed('login', libname);    return;	end;
    disp('login successful');

    % request the tracks available to this user (this is optional)
    numTracks = calllib(libname,'numberOfTracks');     if (numTracks == 0) whenFailed('numberOfTracks', libname);    return;	end;
    disp([num2str(numTracks) ' track(s):']);
    for i=1:numTracks
        trackname = calllib(libname,'trackName',i-1); % indexing from 0
        if (isempty(trackname))         whenFailed('trackName', libname);    return;
        else                            disp([' ' num2str(i) ': ' trackname]);      end;
    end;

    % set the track
    result = calllib(libname,'setTrack',track);
    if (result == 0) whenFailed('setTrack', libname);    return;	end;
    disp(['track set to "' track '"']);

    % obtain number of problems in the track
    numProblems = calllib(libname,'numberOfProblems');
    if (numProblems == 0) whenFailed('numberOfProblems', libname);     return;	end;
    disp(['The track consists of ' num2str(numProblems) ' problems.']);


    % For demonstration purposes we optimize only the first problem in the track.
    for problemID=0:(numProblems-1)
    result = calllib(libname,'setProblem',problemID);
    if (result == 0) whenFailed('setProblem', libname);    return;	end;
    disp(['Problem ID set to ' num2str(problemID)]);

    % obtain problem properties
    dim = calllib(libname,'dimension');
    if (dim == 0) whenFailed('dimension', libname);    return;	end;
    obj = calllib(libname,'numberOfObjectives');
    if (obj == 0) whenFailed('numberOfObjectives', libname);    return;	end;
    bud = calllib(libname,'budget');
    if (bud == 0) whenFailed('budget', libname);    return;	end;
    evals = calllib(libname,'evaluations');
    if (evals < 0) whenFailed('evaluations', libname);    return;	end;
    disp(['problem dimension: ' num2str(dim)]);
    disp(['number of objectives: ' num2str(obj)]);
    disp(['problem budget: ' num2str(bud)]);
    disp(['number of already used up evaluations: ' num2str(evals)]);

    % allocate memory for a search point and function value
    point = zeros(1,dim);
    value = zeros(1,obj);

    % run the optimization loop
    rand('state',123);
    for e=1:bud
        if (e <= evals)
            % If evals > 0 then we have already optimized this problem to some point.
            % Maybe the optimizer crashed or was interrupted.
            %
            % This code demonstrates a primitive recovery approach, namely to replay
            % the history as if it were the actual black box queries. In this example
            % this affects only "bestValue" since random search does not have any
            % state variables.
            % As soon as e >= evals we switch over to black box queries.
            value = libpointer('doublePtr',zeros(1,obj));
            point = libpointer('doublePtr',zeros(1,dim));
        	[result, point, value] = calllib(libname,'history',e-1,point,value);
        	if (result == 0) whenFailed('history', libname);    end;
        else
            % define a search point, here uniformly at random
            for d=1:dim
                point(d) = rand();
            end;
            % query the black box
            value = libpointer('doublePtr',zeros(1,obj));
            [result, TMP, value] = calllib(libname,'evaluate',point,value);
            if (result == 0) whenFailed('evaluate', libname);    return;	end;
        end;

        % In any real algorithm "point" and "value" would update the internals state.
        sx = sprintf('%d ', point);
        sf = sprintf('%d', value);
        disp(['[' num2str(e) '] f(' sx ') = ' sf]);
    end;

    % check that we are indeed done
    evals = calllib(libname,'evaluations');
    if (evals == 0) whenFailed('evaluations', libname);    return;	end;
    if (evals == bud)   disp('optimization finished.');
    else                disp('something went wrong: number of evaluations does not coincide with budget :(');   end;

    end;

    unloadlibrary(libname);

function whenFailed(funcname, libname)
    disp([funcname '() failed: ' calllib(libname,'errorMessage')]);    
    unloadlibrary(libname);
    