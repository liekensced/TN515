
%% -----------------------------------------------------------
% Parameters to set (link budget and environment)
% -----------------------------------------------------------
Pt_dBm = 20;           % Transmit power at antenna input [dBm]
Gt_dBi = 8.71;            % Tx antenna gain toward Rx [dBi]
Gr_dBi = 8.71;            % Rx antenna gain toward Tx [dBi]
f_Hz   = 4.8e9;        % Carrier frequency [Hz] (example: 915 MHz)
ht     = 1.5;          % Tx antenna height above ground [m]
hr     = 1.5;          % Rx antenna height above ground [m]
pol    = 'H';          % Polarization for reflection: 'V' (vertical/TM) or 'H' (horizontal/TE)


% PU synthetic sports surface, low-loss dielectric.
eps_r = 4.5;
sigma = 1e-4;

% Optional extra losses (dB): polarization mismatch, atmospheric, feeder, etc.
ap_dB = 0;  % polarization mismatch loss (Friis-only curve)
a_dB  = 0;  % other lumped losses



%% -----------------------------------------------------------
% Load data (your original code)
% -----------------------------------------------------------
opts = detectImportOptions('measurements.csv');
distVar = opts.VariableNames{1};
rxVar   = opts.VariableNames{2};
opts    = setvartype(opts, {distVar, rxVar}, 'double');
T       = readtable('measurements.csv', opts);
valid   = ~isnan(T.(distVar)) & ~isnan(T.(rxVar));
T       = T(valid, :);

d = T.(distVar)(:);            % distances [m]
Rx_meas_dBm = T.(rxVar)(:);    % measured received power [dBm]

%% -----------------------------------------------------------
% Constants and helper values
% -----------------------------------------------------------
c   = 299792458;                % speed of light [m/s]
lambda = c / f_Hz;              % wavelength [m]
k   = 2*pi/lambda;              % wavenumber [rad/m]

% Convert gains to linear
Gt = 10^(Gt_dBi/10);
Gr = 10^(Gr_dBi/10);

% Critical distance for two-ray cross-over (NS-2 classic)
d_c = 4*pi*ht*hr/lambda;  % [m]

%% -----------------------------------------------------------
% Theoretical models
% -----------------------------------------------------------

% 1) Free-space (Friis) received power in dBm
Lfs_dB = 20*log10(4*pi*d/lambda);  % FSPL
Pr_fs_dBm = Pt_dBm + Gt_dBi + Gr_dBi - Lfs_dB - ap_dB - a_dB;

% 2) Full two-ray with Fresnel ground reflection (complex Γ(θ))
%    Geometry: direct path d1, reflected path d2; grazing angle θ from ground
d1 = sqrt(d.^2 + (ht - hr).^2);
d2 = sqrt(d.^2 + (ht + hr).^2);
theta_graz = atan2(ht+hr, d);                       % angle from ground [rad]
sin_t = sin(theta_graz);
cos_t = cos(theta_graz);

% Complex relative permittivity of ground: εg = εr - j σ/(ωε0)
eps0 = 8.854187817e-12;           % F/m
omega = 2*pi*f_Hz;
eps_g = eps_r - 1i*sigma/(omega*eps0);

% Fresnel reflection coefficient Γ(θ) using grazing-angle form:
% Γ(θ) = (sinθ - X)/(sinθ + X), with:
%   X_h = sqrt(εg - cos^2θ)   (horizontal/TE)
%   X_v = X_h / εg            (vertical/TM)
Xh = sqrt(eps_g - cos_t.^2);
Xv = Xh ./ eps_g;
if upper(pol)=='H'
    X = Xh;
else
    X = Xv;
end
Gamma = (sin_t - X) ./ (sin_t + X);

% Received power via vector sum of rays (Friis factor λ/(4π))
fac = (lambda/(4*pi)).^2;
sumField = exp(-1i*k.*d1)./d1 + Gamma .* exp(-1i*k.*d2)./d2;
Pr_2ray_full_W = (10^(Pt_dBm/10)/1000) * Gt * Gr * fac .* (abs(sumField).^2);
Pr_2ray_full_dBm = 10*log10(Pr_2ray_full_W*1000);

% 3) Asymptotic two-ray (far region, |Γ|≈1, phase ~ destructive)
%    Pr ≈ Pt*Gt*Gr * (ht^2*hr^2 / d^4)   (in linear power)
Pr_2ray_asym_W = (10^(Pt_dBm/10)/1000) * Gt * Gr .* (ht^2 .* hr^2) ./ (d.^4);
Pr_2ray_asym_dBm = 10*log10(Pr_2ray_asym_W*1000);

%% -----------------------------------------------------------
% Simple calibration offsets (optional)
% Models rarely match absolute level exactly; align each model with data by
% adding a single best-fit offset in dB (least-squares)
off_fs   = mean(Rx_meas_dBm - Pr_fs_dBm);
off_2rf  = mean(Rx_meas_dBm - Pr_2ray_full_dBm);
off_2ra  = mean(Rx_meas_dBm - Pr_2ray_asym_dBm);

Pr_fs_dBm_aligned      = Pr_fs_dBm      + off_fs;
Pr_2ray_full_dBm_align = Pr_2ray_full_dBm + off_2rf;
Pr_2ray_asym_dBm_align = Pr_2ray_asym_dBm + off_2ra;

%% -----------------------------------------------------------
% Validation metrics
% -----------------------------------------------------------
rmse = @(x) sqrt(mean((Rx_meas_dBm - x).^2));
RMSE_fs   = rmse(Pr_fs_dBm_aligned);
RMSE_2rf  = rmse(Pr_2ray_full_dBm_align);
RMSE_2ra  = rmse(Pr_2ray_asym_dBm_align);

fprintf('--- Validation ---\n');
fprintf('Critical distance d_c = %.1f m\n', d_c);
fprintf('Alignment offsets: Friis = %.2f dB | Two-ray full = %.2f dB | Two-ray asym = %.2f dB\n', ...
    off_fs, off_2rf, off_2ra);
fprintf('RMSE (aligned): Free-space = %.2f dB | Two-ray (full) = %.2f dB | Two-ray (asym) = %.2f dB\n', ...
    RMSE_fs, RMSE_2rf, RMSE_2ra);

%% -----------------------------------------------------------
% Plots: measurement + models
% -----------------------------------------------------------
figure('Name','Rx vs Distance: measurement +     models');
plot(d, Rx_meas_dBm, 'k-', 'LineWidth', 1.5, 'DisplayName','Measurement'); hold on;
plot(d, Pr_fs_dBm_aligned,      '-','LineWidth',2, 'DisplayName',sprintf('Friis (aligned, +%.1f dB)',off_fs));
plot(d, Pr_2ray_full_dBm_align, '-','LineWidth',2, 'DisplayName',sprintf('Two-ray full (aligned, +%.1f dB)',off_2rf));
plot(d, Pr_2ray_asym_dBm_align, '--','LineWidth',2, 'DisplayName',sprintf('Two-ray asymptotic (aligned, +%.1f dB)',off_2ra));

xlabel('Distance (m)');
ylabel('Rx (dBm)');
grid on;
title('Received power vs distance — measurement and theoretical models');
xline(d_c, ':r', 'd_c (two-ray cross-over)', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
legend('Location','best');
hold off;
saveas(gcf, 'Rx_vs_Distance.png');

% Optional: log-log path-loss slopes visualization
figure('Name','Path-loss trend (log-dB)');
loglog(d, Pt_dBm + Gt_dBi + Gr_dBi - Rx_meas_dBm, 'k-', 'DisplayName','Measured path loss'); hold on;
loglog(d, Pt_dBm + Gt_dBi + Gr_dBi - Pr_fs_dBm_aligned, '-', 'LineWidth',2, 'DisplayName','Friis aligned');
loglog(d, Pt_dBm + Gt_dBi + Gr_dBi - Pr_2ray_full_dBm_align, '-', 'LineWidth',2, 'DisplayName','Two-ray full aligned');
loglog(d, Pt_dBm + Gt_dBi + Gr_dBi - Pr_2ray_asym_dBm_align, '--', 'LineWidth',2, 'DisplayName','Two-ray asym aligned');
grid on;
xlabel('Distance (m) [log]');
ylabel('Path loss (dB)');
title('Path-loss vs distance (log-dB view)');
legend;
saveas(gcf, 'Path_Loss_LogLog.png');

% Free-space only comparison
figure('Name','Free-space Path Loss');
PL_meas = Pt_dBm + Gt_dBi + Gr_dBi - Rx_meas_dBm;
PL_fs = Pt_dBm + Gt_dBi + Gr_dBi - Pr_fs_dBm_aligned;
plot(d, PL_meas, 'k-o', 'LineWidth', 1.5, 'DisplayName','Measured'); hold on;
plot(d, PL_fs, '-', 'LineWidth', 2, 'DisplayName',sprintf('Friis (aligned, +%.1f dB)',off_fs));
xlabel('Distance (m)');
ylabel('Path Loss (dB)');
grid on;
title('Free-space Path Loss: Measurement vs Friis Model');
legend('Location','best');
hold off;
saveas(gcf, 'Free_Space_Path_Loss.png');
