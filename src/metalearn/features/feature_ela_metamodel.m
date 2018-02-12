function ft = feature_ela_metamodel(X, y)
% ft = FEATURE_ELA_METAMODEL(X, y) returns ELA metamodel features
% for dataset [X, y].
%
% Linear and quadratic regression models with or without interactions 
% are fitted to the initial dataset [X, y]. The adjusted coefficient of 
% determination R^2 is returned in each case as an indicator for model 
% accuracy. Features reflecting the size relations of the model 
% coefficients are extracted.
%
% Features:
%   lin_simple.adj_r2          - adjusted R^2 of simple linear model
%   lin_simple.intercept       - intercept of simple linear model
%   lin_simple.coef.min        - minimal absolute value of linear model
%                                coefficients
%   lin_simple.coef.max        - maximal absolute value of linear model
%                                coefficients
%   lin_simple.coef.max_by_min - ratio of max. and min. abs. values of 
%                                linear model coefficients
%   lin_w_interact.adj_r2      - adjusted R^2 of linear model with 
%                                interactions
%   quad_simple.adj_r2         - adjusted R^2 of pure quadratic model 
%   quad_simple.cond           - ratio of max. and min. abs. values of 
%                                quadratic coefficients
%   quad_w_interact.adj_r2     - adjusted R^2 of quadratic model with 
%                                interactions

  if nargin < 2
    help feature_ela_metamodel
    if nargout > 0
      ft = struct();
    end
    return
  end
  
  dim = size(X, 2);
  
  % simple linear model
  lm = fitlm(X, y, 'linear');
  lm_coeff = lm.Coefficients.Estimate;
  ft.lin_simple.adj_r2 = lm.Rsquared.Adjusted;
  ft.lin_simple.intercept = lm_coeff(1);
  ft.lin_simple.coef.min = min(abs(lm_coeff(2:end)));
  ft.lin_simple.coef.max = max(abs(lm_coeff(2:end)));
  ft.lin_simple.coef.max_by_min = max(abs(lm_coeff(2:end))) / min(abs(lm_coeff(2:end)));
  
  % interactions
  im = fitlm(X, y, 'interactions');
  ft.lin_w_interact.adj_r2 = im.Rsquared.Adjusted;
  
  % pure quadratic model
  pqm = fitlm(X, y, 'purequadratic');
  quad_coeff = pqm.Coefficients.Estimate(end-dim+1:end);
  ft.quad_simple.adj_r2 = pqm.Rsquared.Adjusted;
  ft.quad_simple.cond = max(abs(quad_coeff)) / min(abs(quad_coeff));
  
  % quadratic interactions
  X2 = [X, X.^2];
  qim = fitlm(X2, y, 'interactions');
  ft.quad_w_interact.adj_r2 = qim.Rsquared.Adjusted;

end