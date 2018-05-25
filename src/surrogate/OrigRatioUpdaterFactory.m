classdef OrigRatioUpdaterFactory
% TODO
% [ ] remove updaterParams and use DTAdaptive_* parameters instead
  methods (Static)
    function obj = createUpdater(ec, surrogateOpts)
      switch lower(surrogateOpts.updaterType)
        % % these are not supported anymore:
        % case 'rmse'
        %   obj = OrigRatioUpdaterRMSE(surrogateOpts.updaterParams);
        % case 'kendall'
        %   obj = OrigRatioUpdaterKendall(surrogateOpts.updaterParams);
        case 'rankdiff'
          obj = OrigRatioUpdaterRankDiff2(ec, surrogateOpts);
        otherwise
          % including surrogateOpts.updaterType == 'constant'
          %
          % this awfull code is due to backward-compatibility O:-)
          if ~(isfield(surrogateOpts, 'updaterParams')) || isempty(surrogateOpts.updaterParams)
            if ~(isfield(surrogateOpts, 'origEvalsRatio'))
              if ~(isfield(surrogateOpts, 'evoControlRestrictedParam'))
                error('There''s not a parameter for ConstantRatioUpdater');
              else
                p = surrogateOpts.evoControlRestrictedParam;
              end
            else
              p = surrogateOpts.origEvalsRatio;
            end
          else
            p = surrogateOpts.updaterParams;
          end
          
          obj = OrigRatioUpdaterConstant(p);
      end
    end
  end
end
