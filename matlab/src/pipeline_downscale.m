function [lst_downscaled, stats, products] = pipeline_downscale(s2_refl, l8_rad, w_atm, aoi)
% PIPELINE_DOWNSCALE End-to-end multi-sensor thermal downscaling pipeline.
% Fuses Landsat-8 TIRS (thermal) with Sentinel-2 MSI (optical/SWIR) using the MHUTS algorithm.
%
% Syntax:
%   lst_downscaled = pipeline_downscale(s2_refl, l8_rad, w_atm)
%   [lst_downscaled, stats, products] = pipeline_downscale(s2_refl, l8_rad, w_atm, aoi)
%
% Inputs:
%   s2_refl - Sentinel-2 L2A surface reflectance (H_s2 x W_s2 x 9) [B2,B3,B4,B5,B6,B7,B8,B11,B12].
%   l8_rad  - Landsat-8 TIRS Band 10 TOA Radiance (H_l8 x W_l8) [W / (m^2 * sr * um)].
%   w_atm   - Total atmospheric water vapor content [g/cm^2] (e.g. 2.8861 in paper).
%   aoi     - (Optional) Bounding box [row_start, row_end, col_start, col_end] for sub-AOI.
%
% Outputs:
%   lst_downscaled - Sharpened LST raster at 10 m spatial resolution [K].
%   stats          - Statistical accuracy assessment metrics (reproducing Table 3).
%   products       - Struct containing intermediate physical products:
%                    .albedo, .ndvi, .pv, .emi_s2, .lst_l8_30m.
%
% References:
%   Grondona et al. (2022), Elsevier RSASE 26, 100702.

if nargin < 4
    aoi = [];
end

% 1. Derive Sentinel-2 Broadband Albedo (D'Urso & Calera, 2006)
[albedo_s2, idx] = albedo_sentinel(s2_refl, 3, true);
ndvi_s2 = idx.ndvi;

% 2. Derive Sentinel-2 Biophysical Predictors (PV & Emissivity)
[emi_s2, pv_s2] = vegetation_emissivity(ndvi_s2, 'liu_zhang');

% 3. Resample Emissivity to Landsat-8 grid for Radiative Transfer LST
[H_l8, W_l8] = size(l8_rad);
emi_l8 = imresize(emi_s2, [H_l8, W_l8], 'bilinear');

% 4. Retrieve Landsat-8 Land Surface Temperature via GSCM
lst_l8 = lst_gscm(l8_rad, emi_l8, w_atm, 'landsat8');

% 5. Execute MHUTS Spatial Sharpening (30 m -> 10 m)
[lst_downscaled, stats, model] = mhuts_downscale(lst_l8, emi_s2, albedo_s2, pv_s2, aoi);

% 6. Compile Intermediate Products
if nargout >= 3
    products.albedo     = albedo_s2;
    products.ndvi       = ndvi_s2;
    products.pv         = pv_s2;
    products.emi_s2     = emi_s2;
    products.lst_l8     = lst_l8;
    products.model      = model;
end

end
