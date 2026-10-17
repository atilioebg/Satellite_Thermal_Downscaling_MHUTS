function [lst, tb] = lst_gscm(rad, emi, w, sensor)
% LST_GSCM Land Surface Temperature retrieval via Generalized Single-Channel Method (GSCM).
%
% Syntax:
%   lst = lst_gscm(rad, emi, w)
%   [lst, tb] = lst_gscm(rad, emi, w, sensor)
%
% Inputs:
%   rad    - 2D matrix of Top-of-Atmosphere (TOA) thermal radiance [W / (m^2 * sr * um)].
%   emi    - 2D matrix of Land Surface Emissivity (LSE) [0, 1].
%   w      - Total atmospheric water vapor content [g/cm^2] (e.g. from ENVI FLAASH or MODIS).
%   sensor - (Optional) Sensor type string: 'landsat8' [Default, TIRS Band 10].
%
% Outputs:
%   lst    - 2D matrix of Land Surface Temperature in Kelvin [K].
%   tb     - 2D matrix of Top-of-Atmosphere Brightness Temperature in Kelvin [K].
%
% References:
%   Grondona et al. (2022), Elsevier RSASE 26, 100702, Eq (4)-(7).
%   Jimenez-Munoz, J. C., et al. (2014). Land surface temperature retrieval methods 
%   from Landsat 8 TIRS data. IEEE Geoscience and Remote Sensing Letters, 11(10), 1840-1843.

if nargin < 4 || isempty(sensor)
    sensor = 'landsat8';
end

rad = double(rad);
emi = double(emi);
w   = double(w);

switch lower(sensor)
    case 'landsat8'
        % Landsat-8 TIRS Band 10 Calibration Constants
        K1 = 774.89;   % W / (m^2 * sr * um)
        K2 = 1321.08;  % K
        Br = 1324.0;   % K

        % Atmospheric correction coefficients for TIRS Band 10 (Jimenez-Munoz et al., 2014)
        psi1 =  0.04019 * (w^2) + 0.02916 * w + 1.01523;
        psi2 = -0.38333 * (w^2) - 1.50294 * w + 0.20321;
        psi3 =  0.00918 * (w^2) + 1.36072 * w - 0.27514;
        
    otherwise
        error('Unsupported sensor. Currently supporting: ''landsat8''.');
end

% Effective Brightness Temperature (Planck inversion)
tb = K2 ./ log((K1 ./ rad) + 1.0);

% Linearization parameters
gamma = (tb.^2) ./ (Br .* rad);
delta = tb - ((tb.^2) ./ Br);

% Generalized Single-Channel Method (GSCM LST)
lst = gamma .* (((psi1 .* rad + psi2) ./ emi) + psi3) + delta;

% Physical plausibility filter
lst(lst < 200 | lst > 380 | isnan(lst)) = NaN;

end
