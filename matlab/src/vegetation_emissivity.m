function [emi, pv] = vegetation_emissivity(ndvi, method)
% VEGETATION_EMISSIVITY Computes Fractional Vegetation Cover (PV) and Land Surface Emissivity (LSE).
%
% Syntax:
%   [emi, pv] = vegetation_emissivity(ndvi)
%   [emi, pv] = vegetation_emissivity(ndvi, method)
%
% Inputs:
%   ndvi   - 2D matrix (H x W) containing Normalized Difference Vegetation Index.
%   method - (Optional) Emissivity retrieval formulation:
%            'liu_zhang' [Default]: Piecewise semi-empirical model from Liu & Zhang (2011)
%            'sobrino': Thresholded model based on Sobrino et al. (2004)
%
% Outputs:
%   emi    - 2D matrix (H x W) of Land Surface Emissivity (LSE).
%   pv     - 2D matrix (H x W) of Fractional Vegetation Cover [0, 1].
%
% References:
%   Grondona et al. (2022), Elsevier RSASE 26, 100702, Eq (8) & Eq (12).
%   Liu, L., & Zhang, Y. (2011). Urban surface emissivity retrieval.
%   Stathopoulou, M., & Cartalis, C. (2009). Downscaling AVHRR LST.

if nargin < 2 || isempty(method)
    method = 'liu_zhang';
end

ndvi = double(ndvi);

% Fractional Vegetation Cover (PV) - Stathopoulou & Cartalis (2009)
ndvi_min = min(ndvi(:));
ndvi_max = max(ndvi(:));
range_val = ndvi_max - ndvi_min;

if range_val < 1e-6
    pv = zeros(size(ndvi));
else
    pv = ((ndvi - ndvi_min) ./ range_val).^2;
end
pv = max(0.0, min(1.0, pv));

% Land Surface Emissivity (LSE)
emi = zeros(size(ndvi));

switch lower(method)
    case 'liu_zhang'
        % Class 1: Dense Vegetation (NDVI >= 0.727)
        mask_veg = (ndvi >= 0.727);
        emi(mask_veg) = 0.990;

        % Class 2: Mixed Soil-Vegetation Pixels (0.157 <= NDVI < 0.727)
        mask_mix = (ndvi >= 0.157) & (ndvi < 0.727);
        emi(mask_mix) = 1.0094 + 0.047 * log(ndvi(mask_mix));

        % Class 3: Bare Soil / Impervious Surface (-0.185 <= NDVI < 0.157)
        mask_soil = (ndvi >= -0.185) & (ndvi < 0.157);
        emi(mask_soil) = 0.970;

        % Class 4: Water Bodies / Snow (NDVI < -0.185)
        mask_water = (ndvi < -0.185);
        emi(mask_water) = 0.995;

    case 'sobrino'
        emi_soil = 0.97;
        emi_veg = 0.99;
        d_emi = 0.005; % Cavity effect / geometric factor
        
        mask_soil = (ndvi < 0.2);
        mask_veg = (ndvi > 0.5);
        mask_mix = (~mask_soil) & (~mask_veg);
        
        emi(mask_soil) = emi_soil;
        emi(mask_veg) = emi_veg;
        emi(mask_mix) = emi_soil + (emi_veg - emi_soil) .* pv(mask_mix) + d_emi;

    otherwise
        error('Unknown emissivity method. Use ''liu_zhang'' or ''sobrino''.');
end

% Physical bounding
emi = max(0.85, min(1.0, emi));

end
