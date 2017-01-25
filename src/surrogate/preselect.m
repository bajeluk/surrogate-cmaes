function [yModel, x, xValid, z] = preselect(n, cmaesState, model, sampleOpts, varargin)
  if (nargin >= 5 && ~isempty(varargin{1}) && isnumeric(varargin{1}))
    preselectFactor = varargin{1};
  else
    preselectFactor = 50;
  end

  [xLarge, xValidLarge, zLarge] = ...
      sampleCmaesNoFitness(cmaesState.sigma, preselectFactor*cmaesState.lambda, cmaesState, sampleOpts);
  [~, yModelLarge] = model.getModelOutput(xValidLarge');

  [yModelSorted, yId] = sort(yModelLarge);
  x      = xLarge(:,yId(1:n));
  xValid = xValidLarge(:,yId(1:n));
  z      = zLarge(:,yId(1:n));
  yModel = yModelLarge(yId(1:n));
end
