% testModelParamsStatistics
%
% loads *.mat files in the directory 'MATFILESDIR' and calculates statistics
% from the regression data loaded from matrices 'kor':
%   columns  ~ different functions
%   rows     ~ different evaluations thresholds
%   3rd dims ~ different dimensions

MATFILESDIR = '~/prg/surrogate-cmaes/exp/experiments/modeltrain_gp_01';
functions = [1 2 3 5 6 8 10 11 12 13 14 20 21];
dims      = [2 5 10];


% load file-names in the directory MATFILESDIR
files = dir([MATFILESDIR filesep '*.mat']);

% prepare output matrices
meansCorr = zeros(length(files), length(functions), length(dims));
mediansCorr = zeros(length(files), length(functions), length(dims));
stdCorr = zeros(length(files), length(functions), length(dims));
mediansPerSettings = zeros(length(files), length(dims));
bestMedians = zeros(length(files), length(dims));
bestMediansId = zeros(length(files), length(dims));

% process each file
for f = 1:length(files)
  load(files(f).name);

  for dim = 1:length(dims)
    % identify non-NaN entries
    valid = ~isnan(kor(:,:,dim));
    % calculate mean, median a std for each function through eval. thresholds
    for j = 1:size(valid,2)
      meansCorr(f,j,dim) = mean(kor(valid(:,j),j,dim));
      mediansCorr(f,j,dim) = median(kor(valid(:,j),j,dim));
      stdCorr(f,j,dim) = std(kor(valid(:,j),j,dim));
    end
  end
end

% draw boxplots for each file/settings and dimension 
% and calculate medians for each settings and dimension through medians
% calculated in the previous step
figure();
for dim = 1:length(dims)
  subplot(length(dims), 1, dim);
  boxplot(mediansCorr(:,:,dim)');
  title(['correlations (medians from progress) per settings in ' num2str(dims(dim)) 'D']);
  mediansPerSettings(:, dim) = median(mediansCorr(:,:,dim), 2);
  mediansPerSettings(isnan(mediansPerSettings(:,dim)), dim) = -Inf;
  
  [bestMedians(:,dim), bestMediansId(:,dim)] = sort(mediansPerSettings(:,dim), 'descend');
end
