classdef ModelFactory
  methods (Static)
    function obj = createModel(str, modelOptions, xMean)
      switch lower(str)
        case 'gp'
          obj = GpModel(modelOptions, xMean);
        case 'fitrgp'
          obj = GprModel(modelOptions, xMean);
        case 'rf'
          obj = RfModel(modelOptions, xMean);
        case 'bbob'
          obj = PreciseModel(modelOptions, xMean);
        otherwise
          warning(['ModelFactory.createModel: ' str ' -- no such model available']);
          obj = [];
      end
    end
  end
end
