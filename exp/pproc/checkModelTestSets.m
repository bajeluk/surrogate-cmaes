function report = checkModelTestSets(filename)
%checkModelTestSets reports the necessary properties of the model testing
% dataset.
%
% report = checkModelTestSets(filename) - report the properties of the
%   dataset saved in 'filename'.
%
% Input:
%   filename - name of the file with dataset or the dataset itself | string
%              or cell-array of data
%
% Output:
%   report - structure with resulting properties as fields (using the
%            notation of ids: function id - f, dimension id - d, instance
%            id - i, general id - id, generation number - g):
%     .emptySet - N setting ids [f, d, i, id] of empty cells | Nx4 double
%     .testSetX - structure with resulting properties of testSetX as
%                 fields:
%        .wrongSize - N ids [f, d, i, id, g], where the size of testSetX
%                     is not corresponding to the generation lambda | Nx5
%                     double
%        .duplicity - N ids [f, d, i, id, g], where the points are not
%                     unique | Nx5 double
%     .testSetY - structure with resulting properties of testSetY as
%                 fields:
%        .wrongSize - N ids [f, d, i, id, g], where the size of testSetY
%                     is not corresponding to the generation lambda | Nx5
%                     double
%
% See Also:
%   modelTestSets, datasetFromInstances

  % input is a file
  if ischar(filename)
    assert(isfile(filename), 'scmaes:checkModelTestSets:wrongFileName', ...
           'File %s does not exist.', filename)
    % load dataset
    S = load(filename);

  % input is a dataset
  elseif iscell(filename)
    S.ds = filename;
  else
    error('scmaes:checkModelTestSets:wrongInput', ...
          'Input must be name of the file or dataset (cell-array).')
  end

  [nFuns, nDims, nInsts, nIds] = size(S.ds);

  % init checking variables
  emptySet = [];
  wrongSizeTestSetX = [];
  wrongSizeTestSetY = [];
  duplicityTestSetX = [];

  % checking loop
  for f = 1:nFuns
    for d = 1:nDims
      for i = 1:nInsts
        for id = 1:nIds

          if isempty(S.ds{f, d, i, id})
            % save empty ids
            emptySet = [emptySet; f, d, i, id];
          else
            % check test set sizes
            for g = 1:numel(S.ds{f, d, i, id}.generations)
              actualLambda = S.ds{f, d, i, id}.cmaesStates{g}.lambda;
              % verify number of points in testSetX
              if actualLambda ~= size(S.ds{f, d, i, id}.testSetX{g}, 1)
                wrongSizeTestSetX = [wrongSizeTestSetX; f, d, i, id, g];
              end
              % verify number of points in testSetY
              if actualLambda ~= size(S.ds{f, d, i, id}.testSetY{g}, 1)
                wrongSizeTestSetY = [wrongSizeTestSetY; f, d, i, id, g];
              end
              % verify uniqueness of points in testSetX
              if size(unique(S.ds{1}.archive.X, 'rows'), 1) < size(S.ds{1}.archive.X, 1)
                duplicityTestSetX = [duplicityTestSetX; f, d, i, id, g];
              end
            end
          end
        end
      end
    end
  end

  % export all the results
  report.emptySet = emptySet;
  report.testSetX.wrongSize = wrongSizeTestSetX;
  report.testSetX.duplicity = duplicityTestSetX;
  report.testSetY.wrongSize = wrongSizeTestSetY;
end