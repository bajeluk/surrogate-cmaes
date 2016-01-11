classdef ECFactory
  methods (Static)
    function obj = createEC(str, surrogateOpts)
      switch lower(str)
        case 'individual'
          obj = IndividualEC();
        case 'generation'
          obj = GenerationEC(surrogateOpts);
        case {'doubletrained', 'restricted'}
          obj = DoubleTrainedEC();
        otherwise
          warning(['ECFactory.createEC: ' str ' -- no such evolution control available']);
          obj = [];
      end
    end
  end
end
