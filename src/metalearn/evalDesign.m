function Y = evalDesign(X, funId, instId)
%EVALDESIGN Evaluate a sample with a function.
%   Y = evalDesign(X, funId, instId) evaluates all column vectors in
%   design matrix X with BBOB function id funId, instance instId.

   if ~fgeneric('EXIST', funId)
     error('EVALDESIGN: funId not a BBOB function id.');
   end

   bbob_dir = '/tmp/bbob';
   
   if ~exist(bbob_dir, 'dir')
     mkdir(bbob_dir);
   end

   fgeneric('initialize', funId, instId, bbob_dir);

   % scales into BBOB search interval
   Y = fgeneric(X);
   fgeneric('finalize');
   
%    try
%      % cleanup bbob files
%      dat_dir = fullfile(bbob_dir, sprintf('data_f%d', funId));
%      dat_files = fullfile(dat_dir, '*');
%      info_files = fullfile(bbob_dir, '*info');
% 
%      delete(dat_files);
%      delete(info_files);
%      rmdir(dat_dir);
%      rmdir(bbob_dir);
%    catch e
%      getReport(e)
%      warning('Could not clean up %s', bbob_dir);
%    end

end

