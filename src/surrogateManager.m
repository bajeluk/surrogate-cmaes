function [fitness_raw, arx, arxvalid, arz, counteval] = surrogateManager(xmean, sigma, lambda, BD, diagD, fitfun_handle, varargin, surrogateOpts)
  if (strcmpi(func2str(surrogateOpts.sampleFcn), 'samplecmaes'))
    [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, varargin, surrogateOpts.sampleOpts);
  end
end