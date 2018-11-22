function ft = feature_ela_metamodel(X, y, ~)
% ft = FEATURE_ELA_METAMODEL(X, y, settings) returns ELA metamodel features
% for dataset [X, y] according to settings.
%
% Linear and quadratic regression models with or without interactions 
% are fitted to the initial dataset [X, y]. The adjusted coefficient of 
% determination R^2 is returned in each case as an indicator for model 
% accuracy. Features reflecting the size relations of the model 
% coefficients are extracted.
%
% settings: - there are no settings, the field is implemented due to
%             compatibility with the rest of features
%
% Features:
%   lin_simple_adj_r2          - adjusted R^2 of simple linear model
%   lin_simple_intercept       - intercept of simple linear model
%   lin_simple_coef_min        - minimal absolute value of linear model
%                                coefficients
%   lin_simple_coef_max        - maximal absolute value of linear model
%                                coefficients
%   lin_simple_coef_max_by_min - ratio of max. and min. abs. values of 
%                                linear model coefficients
%   lin_w_interact_adj_r2      - adjusted R^2 of linear model with 
%                                interactions
%   quad_simple_adj_r2         - adjusted R^2 of pure quadratic model 
%   quad_simple_cond           - ratio of max. and min. abs. values of 
%                                quadratic coefficients
%   quad_w_interact_adj_r2     - adjusted R^2 of quadratic model with 
%                                interactions

  if nargin < 2
    help feature_ela_metamodel
    if nargout > 0
      ft = struct();
    end
    return
  end
  
  dim = size(X, 2);
  
  % y-values not available
  if all(isnan(y))
    ft.lin_simple_adj_r2 = NaN;
    ft.lin_simple_intercept = NaN;
    ft.lin_simple_coef_min = NaN;
    ft.lin_simple_coef_max = NaN;
    ft.lin_simple_coef_max_by_min = NaN;
    ft.lin_w_interact_adj_r2 = NaN;
    ft.quad_simple_adj_r2 = NaN;
    ft.quad_simple_cond = NaN;
    ft.quad_w_interact_adj_r2 = NaN;
    return
  end
  
  % simple linear model
  lm = fitlm(X, y, 'linear');
  lm_coeff = lm.Coefficients.Estimate;
  ft.lin_simple_adj_r2 = lm.Rsquared.Adjusted;
  ft.lin_simple_intercept = lm_coeff(1);
  ft.lin_simple_coef_min = min(abs(lm_coeff(2:end)));
  ft.lin_simple_coef_max = max(abs(lm_coeff(2:end)));
  ft.lin_simple_coef_max_by_min = max(abs(lm_coeff(2:end))) / min(abs(lm_coeff(2:end)));
  
  % interactions
  im = fitlm(X, y, 'interactions');
  ft.lin_w_interact_adj_r2 = im.Rsquared.Adjusted;
  
  % pure quadratic model
  pqm = fitlm(X, y, 'purequadratic');
  quad_coeff = pqm.Coefficients.Estimate(end-dim+1:end);
  ft.quad_simple_adj_r2 = pqm.Rsquared.Adjusted;
  ft.quad_simple_cond = max(abs(quad_coeff)) / min(abs(quad_coeff));
  
  % quadratic interactions
  X2 = [X, X.^2];
  qim = fitlm(X2, y, 'interactions');
  ft.quad_w_interact_adj_r2 = qim.Rsquared.Adjusted;

  % ensure features to be non-empty in case of empty input
  if isempty(X) || isempty(y)
    ft = repStructVal(ft, @isempty, NaN, 'test');
  end
  
end