classdef ECFactory
  methods (Static)
    function obj = createEC(surrogateOpts)
      switch lower(surrogateOpts.evoControl)
        case 'individual'
          obj = IndividualEC();
        case 'generation'
          obj = GenerationEC(surrogateOpts);
        case {'doubletrained', 'restricted'}
          obj = DoubleTrainedEC();
        case 'none'
          obj = NoneEC();
        otherwise
          warning(['ECFactory.createEC: ', surrogateOpts.evoControl, ' -- no such evolution control available']);
          obj = [];
      end
    end
  end
end
