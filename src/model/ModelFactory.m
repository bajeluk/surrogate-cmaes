classdef ModelFactory
  methods (Static)
    function obj = createModel(str, modelOptions, xMean, oldModel)
      switch lower(str)
        case 'gp'
          obj = GpModel(modelOptions, xMean);
        case 'fitrgp'
          obj = GprModel(modelOptions, xMean);
        case 'rf'
          obj = RfModel(modelOptions, xMean);
        case 'bbob'
          obj = PreciseModel(modelOptions, xMean);
        case 'modelpool'
          % use the supplied 'oldModel' if exists
          if (nargin > 3 && ~isempty(oldModel) && isa(oldModel, 'ModelPool'))
            obj = oldModel;
          else
            obj = ModelPool(modelOptions, xMean);
          end
        otherwise
          warning(['ModelFactory.createModel: ' str ' -- no such model available']);
          obj = [];
      end
    end
  end
end
