function [nF_opt, vX_opt] = bobyqa(sObjFunName, vX_opt, vX_l, vX_u, svOptions)
% The BOBYQA algorithm for bound constrained optimization without
% derivatives by M.J.D. Powell
% 
% This is a sligtly modified version by Lukas Bajer, 2017 (the modifications
% are licensed under CC license). For the original version, see below:
% 
% ==== License ====
% 
% Copyright (c) [2014] [Karlsruhe Institute of Technology
%                       Institute of Engineering Mechanics]
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following condition:
% 
%   * The above copyright notice and this permission notice shall be
%     included in all copies or substantial portions of the Software.
% 
% ==== Preparation ====
% 
% This file uses the dlib C++ implementation of BOBYQA.
% To use the algorithm, the mex file has to be compiled for the specific
% archtecture of your computer. Please refer to the official MATLAB help
% for further details.
% 
% You need three files to run the algorithm:
% 
%   * This file, which executes the mex file, defines default parameters,
%     updates the screen during the iterations and enables the C++
%     algorithm to evaluate MATLAB objective functions.
%   * The file "bobyqa_alg.cpp" which contains the C++ source code defining
%     a gateway between MATLAB and the dlib library containing the
%     algorithm.
%   * The dlib library containing the C++ source files with the BOBYQA
%     algorithm. Successful compilation requires a copy of the (sub)folder
%     "/dlib" in the same directory as the other two files.
%     You can download the dlib library at http://dlib.net
% 
% The function will display the necessary commands for compilation in case
% of an error.
% 
% 
% ==== Execution ====
% 
% This file sets defaults for the BOBYQA parameters (if not specified) and
% executes the mex function. The mex function itself calls this file to
% evaluate the MATLAB objective function during the iterations.
% 
% Input parameters:
% 
%   * sObjFunName   string/function_handle   name of the objective function
% 
%   * vX_opt        vector                   initial optimization vector
% 
%   * vX_l          vector (optional)        lower bounds for optimization
%                                               default: -1e100 (no bound)
% 
%   * vX_u          vector (optional)        upper bounds for optimization
%                                               default:  1e100 (no bound)
% 
%   * svOpts        structure                options for display/algorithm
%           .display         display mode:
%                              'none': no output during iterations
%                              'iter': (default) output after every step
%           .npt             number of points for quadratic approximation
%                              default: 2*n + 1
%           .rho_beg         initial trust region radius
%                              default: 10
%           .rho_beg         final trust region radius
%                              default: 1e-6
%           .maxFunEval      maximum number of function evaluations
%                              default: 1000
% 
% Output parameters:
% 
%   * nF_opt        scalar   value of objective function for final step
% 
%   * vX_opt        vector   final optimization vector
% 
% If the optional parameters are not specified, the default values are set.
% If only one of the bounds (vX_l or vX_u) is specified, the default is set
% for both bounds!
% 
% The file uses a persistent variable to save previous results for display
% purposes. If a critical error occurs during the evaluation and the
% function is not exited properly it could be possible, that the persistent
% variable is not deleted properly. In this case make sure the input
% variable "sObjFunName" is a function handle in the next iteration!
%
%
% ==== Example ====
% 
% To test the algorithm you can run the following MATLAB code as example:
% 
%     fprintf('\n\n>>> First run:\n\n');
%     fhTestfun = @(vX) norm(vX-[3;5;1;7]);
% vX_opt = [-4;5;99;3];
% vX_l = -1e100*ones(4,1);
% vX_u =  1e100*ones(4,1);
% svOptions = struct('display', 'iter', 'npt', 9, 'rho_beg', 10, ...
%                    'rho_end', 1e-6, 'maxFunEval', 100);
% [nF_opt, vX_opt] = bobyqa(fhTestfun, vX_opt, vX_l, vX_u, svOptions)
%     fprintf('\n\n>>> Second run:\n\n');
%     svOptions = struct('display', 'none', 'npt', 9, 'rho_beg', 10, ...
%                        'rho_end', 1e-6, 'maxFunEval', 1000);
%     [nF_opt, vX_opt] = bobyqa(fhTestfun, vX_opt, vX_l, vX_u, svOptions)

narginchk(2, 5);
nargoutchk(0, 2);

% Persistent variable for output during evaluation
persistent svOpts;

try
    if isempty(svOpts) || isa(sObjFunName, 'function_handle') % => first call to function -> initialization
        %% Call BOBYQA_ALG mex function
        nN = numel(vX_opt);
        
        % Set default options if necessary
        svOpts = struct('display', 'iter', ...
                        'npt', 2*nN+1, ...
                        'rho_beg', 10, ...
                        'rho_end', 1e-6, ...
                        'maxFunEval', 1000);
        if nargin == 5
            cOpts = {'display', 'npt', 'rho_beg', 'rho_end', 'maxFunEval'};
            for i = 1:numel(cOpts)
                if isfield(svOptions, cOpts{i})
                    svOpts.(cOpts{i}) = svOptions.(cOpts{i});
                end
            end
        end
        
        % Hold display variables
        svOpts.nFunEval = 0;
        
        if ~exist('vX_l','var') || ~exist('vX_u','var') || isempty(vX_l) || isempty(vX_u)
            vX_l = -1e100*ones(nN,1);
            vX_u =  1e100*ones(nN,1);
        end
        
        if isstr(sObjFunName)
            if exist(sObjFunName, 'file')
                sObjFunName = str2func(sObjFunName);
            else
                error('bobyqa:UnknownFunction', 'The function %s is not a file', sObjFunName);
            end
        end
        
        % Call mex function
        try
            [nF_opt, vX_opt] = dlib_bobyqa(sObjFunName, vX_opt, svOpts.npt, vX_l, vX_u, svOpts.rho_beg, svOpts.rho_end, svOpts.maxFunEval);
        catch ME
            fprintf(['\n\n', 'An error occured trying to evaluate the BOBYQA_ALG mex function.', '\n', ...
                     'If the error perisits try to recompile the mex function for your system.', '\n', ...
                     'Navigate to', '\n\n' ...
                     '    ', strrep(mfilename('fullpath'), '\', '\\'), '\n\n', ...
                     'and use the command', '\n\n', ...
                     '    mex(strcat(''-I"'',pwd,''"''), ''bobyqa_alg.cpp'')', '\n\n']);
            rethrow(ME);
        end
        
    else
        %% Evaluate objective function if called by BOBYQA_ALG
        if nargout < 2 && nargin == 2
            % fhObjFun = str2func(sObjFunName);
            nF_opt = fhObjFun(vX_opt);
            status_display()
            return;
        end
    end
catch ME
    % Clear persistent variable in case of unsuccessful previous run
    clear svOpts;
    rethrow(ME);
end

% Clean up
clear svOpts;

%% Subfunctions
function status_display()
% STATUS_DISPLAY Displays the status in the Matlab command window depending
% on the 'display' option set in the options.

    switch svOpts.display
        case 'none'
        case 'iter'
            if mod(svOpts.nFunEval, 30) == 0
                fprintf(['\n', 'FunEval        ObjFunVal     Norm of step    Rel norm step', '\n']);
            end
            svOpts.nFunEval = svOpts.nFunEval + 1;
            if isfield(svOpts, 'vX_last')
                sNormStep = sprintf('%12.8e', norm(vX_opt - svOpts.vX_last));
                sRelStep  = sprintf('%12.8e', norm(vX_opt - svOpts.vX_last)/svOpts.rho_end);
            else
                sNormStep = '            ';
                sRelStep  = '            ';
            end
            fprintf([sprintf('%7.1u', svOpts.nFunEval), '   ', ...
                     sprintf('%12.8e', nF_opt), '   ', ...
                     sNormStep, '   ', ...
                     sRelStep, ...
                     '\n']);
            svOpts.vX_last = vX_opt;
        otherwise
            fprintf('Unknown display setting. Display set to ''none''.');
            svOpts.display = 'none';
    end
end

end
