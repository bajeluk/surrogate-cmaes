function ft = feature_linear_model(X, y, settings)
% ft = FEATURE_LINEAR_MODEL(X, y, settings) returns cell mapping linear
% model features for dataset [X, y] according to settings.
%
% The cell mapping features discretizes the continuous input space
% utilizing a pre-defined number of blocks (cells) per dimension. Linear 
% model features fit a linear model within each cell of the feature object 
% and aggregate the coefficient vectors across all the cells. 
% (Kerschke et al., 2017)
%
% settings:
%   blocks   - number of cell blocks per dimension
%   lb       - lower bounds of the input space
%   ub       - upper bounds of the input space
%
% Features:
%   lm_avg_length_reg  - length of the average coefficient vector
%   lm_avg_length_norm - length of the average normalized coefficient 
%                        vector
%   lm_length_mean     - mean of the lengths of the coefficient vectors
%   lm_length_std      - standard deviation of the lengths of the 
%                        coefficient vectors
%   lm_corr_reg        - correlation of the coefficient vectors
%   lm_corr_norm       - correlation of the normalized coefficient vectors
%   lm_ratio_mean      - mean of the ratio of maximum and minimum of the
%                        (non-intercept) coefficients
%   lm_ratio_std       - standard deviation of the ratio of maximum and 
%                        minimum of the (non-intercept) coefficients
%   lm_std_radio_reg   - ratio of maximum and minimum of the standard 
%                        deviation of the (non-intercept) coefficients
%   lm_std_radio_norm  - ratio of maximum and minimum of the standard 
%                        deviation of the (non-intercept) normalized 
%                        coefficients
%   lm_std_mean_reg    - mean of the standard deviation of the 
%                        (non-intercept) coefficients
%   lm_std_mean_norm   - mean of the standard deviation of the 
%                        (non-intercept) normalized coefficients 

  if nargin < 3
    if nargin < 2
      help feature_linear_model
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end

  % parse settings
  lb = defopts(settings, 'lb', min(X));
  ub = defopts(settings, 'ub', max(X));
  blocks = defopts(settings, 'blocks', 2);
  
  % create grid of cells
  cmg = CMGrid(X, y, lb, ub, blocks);
  
  % fit linear model in each cell
  lm = cmg.fitPolyModel('linear');
  % get coefficients without intercept
  lm_coeff = NaN(cmg.nCells, cmg.dim);
  for c = 1:cmg.nCells
    if ~isempty(lm{c}.Coefficients)
      lm_coeff(c, :) = lm{c}.Coefficients.Estimate(2:end);
    end
  end
  % remove NaN rows (cells with not enough data to train the model)
  lm_coeff(all(isnan(lm_coeff), 2), :) = [];
  nCoeff = size(lm_coeff, 1);
  
  sumCells = prod(cmg.blocks);
  % warn in case of empty cells
  if nCoeff < sumCells
    warning(['%d out of %d cells (%0.2f%%) contains too few points to construct a linear model.', ...
             'This may make linear model features useless.'], ...
            sumCells - nCoeff, sumCells, (sumCells - nCoeff)/sumCells * 100)
  end
  
  % in case of too few points in all cells return NaN
  if nCoeff == 0
    lm_coeff = NaN(1, cmg.dim);
    nCoeff = 1;
  end
  
  % ratio of biggest to smallest absolute coefficient (per cell)
  lm_coeff_ratio = max(abs(lm_coeff), [], 2) ./ min(abs(lm_coeff), [], 2);
  % normalized vectors of coefficients
  lm_coeff_norm = lm_coeff ./ norm(lm_coeff);
  % length of each single coefficient vector
  lm_coeff_length = sqrt(sum(lm_coeff.^2, 2));
  % standard deviation of coefficients (per feature)
  lm_coeff_std = std(lm_coeff, 1);
  % standard deviation of normalized coefficients (per feature)
  lm_coeff_norm_std = std(lm_coeff_norm, 1);
  % correlation between coefficients
  lm_coeff_corr = corr(lm_coeff') + diag(NaN(nCoeff, 1));
  lm_coeff_corr = sum(sum(lm_coeff_corr, 'omitnan')) / (nCoeff*(nCoeff-1));
  % correlation between normalized coefficients
  lm_coeff_norm_corr = corr(lm_coeff_norm') + diag(NaN(nCoeff, 1));
  lm_coeff_norm_corr = sum(sum(lm_coeff_norm_corr, 'omitnan')) / (nCoeff*(nCoeff-1));
  
  % calculate features
  ft.lm_avg_length_reg  = sqrt(sum(mean(lm_coeff, 2).^2));
  ft.lm_avg_length_norm = sqrt(sum(mean(lm_coeff_norm, 2).^2));
  ft.lm_length_mean = mean(lm_coeff_length);
  ft.lm_length_std  = std(lm_coeff_length);
  ft.lm_corr_reg  = lm_coeff_corr;
  ft.lm_corr_norm = lm_coeff_norm_corr;
  ft.lm_ratio_mean = mean(lm_coeff_ratio);
  ft.lm_ratio_std  = std(lm_coeff_ratio);
  ft.lm_std_radio_reg  = max(lm_coeff_std) / min(lm_coeff_std);
  ft.lm_std_radio_norm = max(lm_coeff_norm_std) / min(lm_coeff_norm_std);
  ft.lm_std_mean_reg  = mean(lm_coeff_std);
  ft.lm_std_mean_norm = mean(lm_coeff_norm_std);
  
end