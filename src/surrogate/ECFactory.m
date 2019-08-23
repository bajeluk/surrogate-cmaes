classdef ECFactory
  methods (Static)
    function obj = createEC(surrogateOpts)
      switch lower(surrogateOpts.evoControl)
        case 'individual'
          obj = IndividualEC();
        case 'generation'
          obj = GenerationEC(surrogateOpts);
        case {'doubletrained', 'doublytrained', 'restricted'}
          obj = DoubleTrainedEC(surrogateOpts);
        case 'maes'
          obj = ModelAssistedEC();
        case 'multitrained'
          obj = MultiTrainedEC(surrogateOpts);
        case 'none'
          obj = NoneEC(surrogateOpts);
        case 'linquad'
          obj = LinQuadEC(surrogateOpts);
        otherwise
          warning(['ECFactory.createEC: ', surrogateOpts.evoControl, ' -- no such evolution control available']);
          obj = [];
      end
    end
  end
end
