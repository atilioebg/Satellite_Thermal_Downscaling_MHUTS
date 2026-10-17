function [stats, stat_table] = evaluate_metrics(y_true, y_pred, label)
% EVALUATE_METRICS Computes validation statistics and formats metrics as in Table 3.
%
% Syntax:
%   stats = evaluate_metrics(y_true, y_pred)
%   [stats, stat_table] = evaluate_metrics(y_true, y_pred, label)
%
% Inputs:
%   y_true - Ground truth or reference coarse LST array [K].
%   y_pred - Predicted / downscaled high-resolution LST array [K].
%   label  - (Optional) String descriptor for the evaluation zone (e.g. 'Landfill AOI').
%
% Outputs:
%   stats      - Struct containing fields: .mean, .std, .rmse, .max, .min.
%   stat_table - Formatted table reproducing the 5 columns of Table 3 in Grondona et al. (2022).
%
% Reference:
%   Grondona et al. (2022), Elsevier RSASE 26, 100702, Table 3.

if nargin < 3 || isempty(label)
    label = 'Evaluation Zone';
end

yt = y_true(:);
yp = y_pred(:);

% Remove NaNs if any
valid = ~isnan(yt) & ~isnan(yp);
yt = yt(valid);
yp = yp(valid);

stats.label = label;
stats.mean  = mean(yp);
stats.std   = std(yp);
stats.rmse  = sqrt(mean((yt - yp).^2));
stats.max   = max(yp);
stats.min   = min(yp);

% Create structured Table (identical to Table 3 in the Elsevier paper)
stat_table = table({label}, stats.mean, stats.std, stats.rmse, stats.max, stats.min, ...
    'VariableNames', {'Scope', 'Mean_K', 'StdDev_K', 'RMSE_K', 'Max_K', 'Min_K'});

end
