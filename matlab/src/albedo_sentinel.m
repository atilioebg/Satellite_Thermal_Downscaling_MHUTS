function [albedo_img, indices] = albedo_sentinel(reflec, method, compute_indices)
% ALBEDO_SENTINEL Computes broadband surface albedo from Sentinel-2 MSI L2A surface reflectance.
%
% Syntax:
%   albedo_img = albedo_sentinel(reflec)
%   albedo_img = albedo_sentinel(reflec, method)
%   [albedo_img, indices] = albedo_sentinel(reflec, method, compute_indices)
%
% Inputs:
%   reflec          - 3D matrix (H x W x 9) containing Sentinel-2 surface reflectance bands:
%                     Band 1: B2  (Blue, 490 nm)
%                     Band 2: B3  (Green, 560 nm)
%                     Band 3: B4  (Red, 665 nm)
%                     Band 4: B5  (Red Edge 1, 705 nm)
%                     Band 5: B6  (Red Edge 2, 740 nm)
%                     Band 6: B7  (Red Edge 3, 783 nm)
%                     Band 7: B8  (NIR, 842 nm)
%                     Band 8: B11 (SWIR 1, 1610 nm)
%                     Band 9: B12 (SWIR 2, 2190 nm)
%   method          - (Optional) Albedo integration method:
%                     1 = Knap et al. (1999)
%                     2 = Liang (2001)
%                     3 = D'Urso & Calera (2006) [Default, calibrated for paper]
%   compute_indices - (Optional) Logical flag (true/false) to return spectral indices
%                     (NDVI, NDBI, NDWI). Default is false.
%
% Outputs:
%   albedo_img      - 2D matrix (H x W) containing broadband surface albedo [0, 1].
%   indices         - (Optional) Struct containing .ndvi, .ndbi, .ndwi.
%
% Reference:
%   Grondona et al. (2022), Elsevier RSASE 26, 100702.
%   D'Urso, G., & Calera, A. (2006). Oper. Meth. to Calc. Crop Water Req.

if nargin < 2 || isempty(method)
    method = 3; % D'Urso and Calera (2006)
end
if nargin < 3
    compute_indices = false;
end

reflec = double(reflec);
if max(reflec(:)) > 1.0
    reflec = reflec / 10000.0; % Scale standard ESA S2 L2A integers to [0, 1]
end

B2  = reflec(:, :, 1);
B3  = reflec(:, :, 2);
B4  = reflec(:, :, 3);
B5  = reflec(:, :, 4);
B6  = reflec(:, :, 5);
B7  = reflec(:, :, 6);
B8  = reflec(:, :, 7);
B11 = reflec(:, :, 8);
B12 = reflec(:, :, 9);

switch method
    case 1
        % Knap et al. (1999)
        albedo_img = 0.726 * B3 - 0.0322 * (B3.^2) - 0.015 * B8 + 0.581 * (B8.^2);
        
    case 2
        % Liang (2001)
        albedo_img = 0.356 * B2 + 0.130 * B4 + 0.373 * B8 + 0.085 * B11 + 0.072 * B12 - 0.0018;
        
    case 3
        % D'Urso & Calera (2006) - 9-band weighted integration (Eq 11 in Grondona et al., 2022)
        albedo_img = 0.1324 * B2 + ...
                     0.1269 * B3 + ...
                     0.1051 * B4 + ...
                     0.0971 * B5 + ...
                     0.0890 * B6 + ...
                     0.0818 * B7 + ...
                     0.0722 * B8 + ...
                     0.0167 * B11 + ...
                     0.0002 * B12;
    otherwise
        error('Invalid method: choose 1 (Knap), 2 (Liang), or 3 (D''Urso).');
end

% Physical bounding
albedo_img = max(0.0, min(1.0, albedo_img));

if compute_indices || nargout > 1
    eps_val = 1e-6;
    indices.ndvi = (B8 - B4) ./ (B8 + B4 + eps_val);
    indices.ndbi = (B11 - B8) ./ (B11 + B8 + eps_val);
    indices.ndwi = (B3 - B8) ./ (B3 + B8 + eps_val);
end

end
