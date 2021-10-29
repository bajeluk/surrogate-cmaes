classdef ModelFactory
  methods (Static)
    function obj = createModel(str, modelOptions, xMean, oldModel)
      switch lower(str)
        case 'gp'
          % Check whether we have defined some parameters on per-dimension
          % basis, i.e. in a cellarray where each cell corresponds to
          % one dimension from 'parameterSets_dimensions' field
          if (isfield(modelOptions, 'parameterSets_dimensions') ...
              && isnumeric(modelOptions.parameterSets_dimensions))
            % Identify the right settings according to current dimension.
            % Dimensions for which we have exact parameterSets are
            % saved in 'modelOptions.parameterSets_dimensions'
            dimensions = modelOptions.parameterSets_dimensions;
            fields = fieldnames(modelOptions);
            dim = size(xMean, 2);
            [~, idDimForParamsPerDim] = min(abs( ...
                modelOptions.parameterSets_dimensions - dim ));
            for fi = 1:length(fields)
              if (iscell(modelOptions.(fields{fi})) ...
                  && length(modelOptions.(fields{fi})) == length(dimensions))
                modelOptions.(fields{fi}) = modelOptions.(fields{fi}){idDimForParamsPerDim};
              end
            end
          end
          obj = GpModel(modelOptions, xMean);
        case 'fitrgp'
          obj = GprModel(modelOptions, xMean);
        case 'rf'
          obj = RfModel(modelOptions, xMean);
        case 'xgb'
          obj = ForestModel(modelOptions, xMean);
        case 'forest'
          obj = ForestModel(modelOptions, xMean);
        case {'bbob', 'precise'}
          obj = PreciseModel(modelOptions, xMean);
        case 'anticorr'
          obj = AntiCorrelatedModel(modelOptions, xMean);
        case 'random'
          obj = RandomModel(modelOptions, xMean);
        case 'reg'
          obj = RegModel(modelOptions, xMean);
        case 'lmm'
          obj = LmmModel(modelOptions, xMean);
        case {'lq', 'hansen'}
          obj = HansenModel(modelOptions, xMean);
        case 'modelpool'
          % use the supplied 'oldModel' if exists
          if (nargin > 3 && ~isempty(oldModel) && isa(oldModel, 'ModelPool'))
            obj = oldModel;
          else
            if (isfield(modelOptions, 'parameterSets_fullfact') ...
                && isnumeric(modelOptions.parameterSets))
              % Identify the right settings according to current dimension.
              % Dimensions for which we have exact parameterSets are
              % saved in 'modelOptions.parameterSets_dimensions'
              dim = size(xMean, 2);
              [~, fullfactIndex] = min(abs( ...
                  modelOptions.parameterSets_dimensions - dim ));
              modelOptions.parameterSets = modelOptions.parameterSets_fullfact( ...
                  modelOptions.parameterSets(fullfactIndex, :) );
            end
            obj = ModelPool(modelOptions, xMean);
          end
        case {'meta', 'metamodel'}
          % use the supplied 'oldModel' if exists
          if (nargin > 3 && ~isempty(oldModel) && isa(oldModel, 'MetaModel'))
            obj = oldModel;
          else
            if (isfield(modelOptions, 'parameterSets_fullfact') ...
                && isnumeric(modelOptions.parameterSets))
              % Identify the right settings according to current dimension.
              % Dimensions for which we have exact parameterSets are
              % saved in 'modelOptions.parameterSets_dimensions'
              dim = size(xMean, 2);
              [~, fullfactIndex] = min(abs( ...
                  modelOptions.parameterSets_dimensions - dim ));
              modelOptions.parameterSets = modelOptions.parameterSets_fullfact( ...
                  modelOptions.parameterSets(fullfactIndex, :) );
            end
            obj = MetaModel(modelOptions, xMean);
          end
        otherwise
          warning(['ModelFactory.createModel: ' str ' -- no such model available']);
          obj = [];
      end
    end
  end
end
