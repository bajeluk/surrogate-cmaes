function [ft, values] = getMetaFeatures(X, y, settings)
% ft = getMetaFeatures(X, y, settings) calculates all implemented
% metafeature groups on dataset [X, y], where X is NxM double matrix of 
% datapoints in the input space and y is Nx1 double vector of objective
% values.
%
% [ft, values] = getMetaFeatures(...) returns also Kx1 vector of
% metafeature values, where K is overall number of metafeatures
%
% Metafeature groups:
%   basic            - features of the initial design
%   cm_angle         - information based on the location of the best and 
%                      worst point of the cell (cell mapping features)
%   cm_convexity     - convexity based on representative observations of 
%                      successive cells (cell mapping features)
%   cm_gradhomo      - cell-wise information on the gradients between each 
%                      point of a cell and its corresponding nearest 
%                      neighbor (cell mapping features)
%   cmaes            - CMA-ES state variable features
%   dispersion       - dispersion among observations within the initial 
%                      design and among a subset of these points
%   ela_distribution - estimation of the density of the initial design’s 
%                      objective values
%   ela_levelset     - discriminant analysis is used to predict groups of
%                      y-values splitted previously by certain treshold
%   ela_metamodel    - regression models are fitted to the initial dataset
%   gcm              - generalized cell mapping: transition probability for
%                      moving from one cell to one of its neighboring cells
%   infocontent      - information content of a continuous landscape
%   linear_model     - statistics of linear model coefficients within each 
%                      cell
%   nearest_better   - comparison of the sets of distances from all 
%                      observations towards their nearest neighbors
%   pca              - principal component analysis on the initial design
%
% Input:
%   X        - values in input space | NxM double matrix
%   y        - objective values | Nx1 double vector
%   settings - structure containing list of feature groups and its 
%              individual settings
%     overall:
%       features  - list of feature groups to calculate | cell-array or 
%                   'all' | {'basic', 'cm_angle', 'cm_convexity', 
%                   'cm_gradhomo', 'dispersion', 'ela_distribution', 
%                   'ela_levelset', 'ela_metamodel', 'gcm', 'infocontent', 
%                   'linear_model', 'nearest_better', 'pca'}
%       lb        - lower bounds of the input space | 1xM double
%       ub        - upper bounds of the input space | 1xM double
%       [feature] - local setting for individual feature group | struct
%       [other]   - fields different from the listed above will be
%                   considered as the settings for all features (local
%                   settings is more important than global)
%     
% Output:
%   ft - structure of metafeatures (first level - groups, second level -
%        metafeatures
%
% Example:
%   X = rand(100, 2);
%   y = rand(100, 1);
%   settings.lb = [0, 0];
%   settings.ub = [1, 1];
%   settings.features = {'cm_convexity', 'cm_gradhomo'};
%   settings.blocks = [5 4];              % settings for all feature groups
%   settings.cm_gradhomo.blocks = [3 4];  % blocks for cm_gradhomo will be
%                                         % different from global setting
%   ft = getMetaFeatures(X, y, settings)
%
%   ft = 
% 
%       cm_convexity: [1x1 struct]
%        cm_gradhomo: [1x1 struct]
%
%   printStructure(ft)
%
%   ft.cm_convexity.concave_soft = 0.636364;
%   ft.cm_convexity.concave_hard = 0.409091;
%   ft.cm_convexity.convex_soft = 0.363636;
%   ft.cm_convexity.convex_hard = 0.318182;
%   ft.cm_gradhomo.grad_mean = 0.326304;
%   ft.cm_gradhomo.grad_std = 0.261862;
%
% See Also:
%   feature_basic
%   feature_cm_angle
%   feature_cm_convexity
%   feature_cm_gradhomo
%   feature_cmaes
%   feature_dispersion
%   feature_ela_distribution
%   feature_ela_levelset
%   feature_ela_metamodel
%   feature_gcm
%   feature_infocontent
%   feature_linear_model
%   feature_nearest_better
%   feature_pca
  
  values = [];
  if nargin < 3
    if nargin < 2
      help getMetaFeatures
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end

  % parse settings
  lb = defopts(settings, 'lb', min(X) - eps);
  ub = defopts(settings, 'ub', max(X) + eps);
  features = defopts(settings, 'features', 'all');
  listFeatures = {'basic', ...
                  'cmaes', ...
                  'cm_angle', ...
                  'cm_convexity', ...
                  'cm_gradhomo', ...
                  'cmaes', ...
                  'dispersion', ...
                  'ela_distribution', ...
                  'ela_levelset', ...
                  'ela_metamodel', ...
                  'gcm', ...
                  'infocontent', ...
                  'linear_model', ...
                  'nearest_better', ...
                  'pca' ...
                 };
  if strcmp(features, 'all')
    features = listFeatures;
  else
    if ~iscell(features)
      features = {features};
    end
    assert(all(cellfun(@(x) any(strcmp(x, listFeatures)), features)), ...
      'Feature names are not correct.')
    assert(numel(features) > 0, 'No features were selected')
  end
  nFeat = numel(features);
  
  % create base structure for individual group settings
  dim = size(X, 2);
  feat_settings_base.lb = myeval(lb);
  feat_settings_base.ub = myeval(ub);
  % add additional overall settings
  commonSettingsList = [{'features', 'lb', 'ub'}, listFeatures];
  set_fields = fieldnames(settings);
  for sf = 1:numel(set_fields)
    if ~any(strcmp(set_fields{sf}, commonSettingsList))
      feat_settings_base.(set_fields{sf}) = settings.(set_fields{sf});
    end
  end
  
  % features loop
  for f = 1:nFeat
    % get overall settings
    feat_settings = feat_settings_base;
    % add feature group specific settings
    group_settings = defopts(settings, features{f}, struct());
    group_fields = fieldnames(group_settings);
    for gf = 1 : numel(group_fields)
      feat_settings.(group_fields{gf}) = group_settings.(group_fields{gf});
    end
    
    % uncomment for debugging:
    % tic
    % fprintf('feature_%s\n', features{f})
    % run feature group calculation
    ft.(features{f}) = eval(['feature_', features{f}, '(X, y, feat_settings)']);
    % toc

    % struct2array unavailable since 2017a
    % see: https://www.mathworks.com/matlabcentral/answers/357399-conversion-of-structure-to-double-in-r2017a
    ft1 = ft.(features{f});
    c = struct2cell(ft1);
    values = [values; [c{:}]'];
  end

end