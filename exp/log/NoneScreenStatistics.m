classdef NoneScreenStatistics < Observer
%SCREENSTATISTICS -- print statistics from DoubleTrainEC on screen without adaptivity
  properties
    verbosity
  end

  methods
    function obj = NoneScreenStatistics(params)
      obj@Observer();
      verbosity = defopts(params, 'verbose', 5);
    end

    function notify(obj, ec, varargin)
      % get the interesting data and process them

      if (~ isfield(ec.surrogateOpts.modelOpts, 'bbob_func'))
        ec.surrogateOpts.fopt = 0;
      end
      if (mod(ec.cmaesState.countiter, 10) == 1)
      %           #####  iter /evals(or:p,b) | Dopt |rmseR | rnkR | rnk2 |rnkVal * | Mo nD nDiR |sigm^2| aErr |smooEr| orRat| aGain|
        fprintf('####### iter /evals | D_fopt. | sigma^2. \n');
      end
      model = '.';
      nTrainData = 0;
      outputValues1 = [...
          ec.cmaesState.countiter, ec.counteval, ...
          ec.stats.fmin - ec.surrogateOpts.fopt];
      outputValues2 = [ ec.cmaesState.sigma^2 ];
      outputValues1(isnan(outputValues1)) = 0.0;
      outputValues2(isnan(outputValues2)) = 0.0;
      %         #####  iter /evals(or,p) | Dopt |rmseR | rnkR | rnk2 |rnkVal * | Mo nD nDiR |sigm^2| aErr |smooEr| orRat| aGain|
      fprintf('=[DTS]= %4d /%5d | %.1e | %.1e \n', outputValues1(:), outputValues2(:) );
    end
  end
end
