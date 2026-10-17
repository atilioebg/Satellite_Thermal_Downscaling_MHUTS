function [lst_downscaled, stats, model] = mhuts_downscale(lst_coarse, emi_fine, albedo_fine, pv_fine, aoi)
% MHUTS_DOWNSCALE Modified High-resolution Urban Thermal Sharpener (MHUTS).
% Downscales thermal infrared Land Surface Temperature (LST) from coarse resolution
% (e.g., Landsat-8 at 30 m) to fine resolution (e.g., Sentinel-2 at 10 m) using
% multivariate cubic polynomial regression against biophysical predictors (PV, Emissivity, Albedo).
%
% Syntax:
%   lst_downscaled = mhuts_downscale(lst_coarse, emi_fine, albedo_fine, pv_fine)
%   [lst_downscaled, stats] = mhuts_downscale(lst_coarse, emi_fine, albedo_fine, pv_fine, aoi)
%   [lst_downscaled, stats, model] = mhuts_downscale(...)
%
% Inputs:
%   lst_coarse  - 2D matrix of coarse-scale LST [K] (e.g. Landsat-8 30 m grid).
%   emi_fine    - 2D matrix of fine-scale Surface Emissivity (Sentinel-2 10 m).
%   albedo_fine - 2D matrix of fine-scale Broadband Albedo (Sentinel-2 10 m).
%   pv_fine     - 2D matrix of fine-scale Fractional Vegetation Cover (Sentinel-2 10 m).
%   aoi         - (Optional) Bounding box [row_start, row_end, col_start, col_end]
%                 to subset the region of interest (e.g. [640, 960, 500, 760] for 
%                 the São Leopoldo Landfill AOI). Default is full extent.
%
% Outputs:
%   lst_downscaled - 2D matrix of sharpened, high-resolution LST [K] at fine scale.
%   stats          - Struct with validation metrics (.mean, .std, .rmse, .max, .min).
%   model          - Trained MATLAB LinearModel object containing the 20 polynomial coefficients.
%
% References:
%   Grondona et al. (2022), Elsevier RSASE 26, 100702, Eq (13) & Table 3.

% 1. Spatial Alignment & Upsampling
[H_fine, W_fine] = size(emi_fine);
[H_coarse, W_coarse] = size(lst_coarse);

if H_coarse ~= H_fine || W_coarse ~= W_fine
    % Compute spatial scaling factor (nominally 3x for 30m -> 10m)
    scale_r = ceil(H_fine / H_coarse);
    scale_c = ceil(W_fine / W_coarse);
    
    % Vectorized Kronecker replication
    lst_upsampled = kron(lst_coarse, ones(scale_r, scale_c));
    
    % Align exact bounding boundaries
    lst_upsampled = lst_upsampled(1:H_fine, 1:W_fine);
else
    lst_upsampled = lst_coarse;
end

% 2. Sub-AOI Extraction (if specified)
if nargin >= 5 && ~isempty(aoi)
    r_start = aoi(1);
    r_end   = aoi(2);
    c_start = aoi(3);
    c_end   = aoi(4);
    
    pv_sub     = pv_fine(r_start:r_end, c_start:c_end);
    emi_sub    = emi_fine(r_start:r_end, c_start:c_end);
    alb_sub    = albedo_fine(r_start:r_end, c_start:c_end);
    lst_up_sub = lst_upsampled(r_start:r_end, c_start:c_end);
    scope_name = 'Sub-AOI Target';
else
    pv_sub     = pv_fine;
    emi_sub    = emi_fine;
    alb_sub    = albedo_fine;
    lst_up_sub = lst_upsampled;
    scope_name = 'Full Extent';
end

% 3. Vectorized Column Reformulation
p = double(pv_sub(:));
e = double(emi_sub(:));
a = double(alb_sub(:));
y = double(lst_up_sub(:));

% Filter invalid / nodata pixels if present
valid = ~isnan(p) & ~isnan(e) & ~isnan(a) & ~isnan(y) & (y > 0);
p = p(valid);
e = e(valid);
a = a(valid);
y = y(valid);

% 4. Least-Squares Multivariate Polynomial Fit ('poly333' -> 20 parameters)
% Eq (13): LST = sum_{i+j+k <= 3} P_{ijk} * (PV)^i * (emi)^j * (albedo)^k
model = fitlm([p, e, a], y, 'poly333');

% 5. Spatial Prediction at Fine Resolution
t_pred = predict(model, [pv_sub(:), emi_sub(:), alb_sub(:)]);
lst_downscaled = reshape(t_pred, size(pv_sub));

% 6. Statistical Validation
if nargout >= 2
    stats = evaluate_metrics(lst_up_sub, lst_downscaled, scope_name);
end

end
