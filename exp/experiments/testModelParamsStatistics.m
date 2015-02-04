% testModelParamsStatistics
%
% loads *.mat files in the directory 'MATFILESDIR' and calculates statistics
% from the regression data loaded from matrices 'kor':
%   columns  ~ different functions
%   rows     ~ different evaluations thresholds
%   3rd dims ~ different dimensions

MATFILESDIR = fullfile('exp', 'experiments', 'modeltrain_gp_01');
functions = [1 2 3 5 6 8 10 11 12 13 14 20 21];
dims      = [2 5 10];


% load file-names in the directory MATFILESDIR
files = dir([MATFILESDIR filesep '*.mat']);

fileIdxs = 1:length(files);
% chosen good GP model settings:
% fileIdxs = [3 4 7 8 11 23 32 39 40 43 44 47 48 63 75 76 79 87 88 95 99]; 

nFiles = length(fileIdxs);

% prepare output matrices
meansCorr = zeros(nFiles, length(functions), length(dims));
mediansCorr = zeros(nFiles, length(functions), length(dims));
stdCorr = zeros(nFiles, length(functions), length(dims));
mediansPerSettings = zeros(nFiles, length(dims));
bestMedians = zeros(nFiles, length(dims));
bestMediansId = zeros(nFiles, length(dims));

% process each file with indeces in 'fileIdxs'
fileNum = 1;
for f = fileIdxs
  load(files(f).name, 'kor');

  for dim = 1:length(dims)
    % identify non-NaN entries
    valid = ~isnan(kor(:,:,dim));
    % calculate mean, median a std for each function through eval. thresholds
    for j = 1:size(valid,2)
      meansCorr(fileNum,j,dim) = mean(kor(valid(:,j),j,dim));
      mediansCorr(fileNum,j,dim) = median(kor(valid(:,j),j,dim));
      stdCorr(fileNum,j,dim) = std(kor(valid(:,j),j,dim));
    end
  end
  fileNum = fileNum + 1;
end

% prepare tick-labels (reduce their frequency for nFiles > 30)
figure();
if (nFiles > 30 && all(fileIdxs == 1:nFiles))
  tickLabels = cell(1,nFiles);
  tickLabels = cellMap(tickLabels, @(x) ' ');
  ticks = [1 5:5:nFiles];
  for i = 1:length(ticks)
    tickLabels{ticks(i)} = num2str(ticks(i));
  end
else
  tickLabels = fileIdxs;
end

% draw boxplots for each file/settings and dimension 
% and calculate medians for each settings and dimension through medians
% calculated in the previous step
for dim = 1:length(dims)
  subplot(length(dims), 1, dim);
  boxplot(mediansCorr(:,:,dim)');
  title(['correlations (medians from progress) per settings in ' num2str(dims(dim)) 'D']);
  set(gca, 'XTick', 1:nFiles);
  set(gca, 'XTickLabel', tickLabels);
  mediansPerSettings(:, dim) = median(mediansCorr(:,:,dim), 2);
  mediansPerSettings(isnan(mediansPerSettings(:,dim)), dim) = -Inf;
  
  [bestMedians(:,dim), bestMediansId(:,dim)] = sort(mediansPerSettings(:,dim), 'descend');
end
