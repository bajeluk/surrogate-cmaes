function hash = modelHash(modelOptions)
% produces short numerical string based on struct modelOptions with
% settings of a model. It uses printStructure() for conversion struct
% into a string
%
% TODO: proper model hash
  if (isempty(modelOptions) || ~isstruct(modelOptions))
    hash = '0';
    return;
  end

  S = printStructure(modelOptions, 'Format', 'field');
  S = double(S);
  % exclude not necessary characters
  S = S(S > 32 & S~= 61) - 32;

  % create hash
  hash = num2str(sum(S.*(1:length(S))));
end
