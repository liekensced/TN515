
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
d_c = 20*pi*ht*hr/lambda/3;  % [m]
d_break = 2*pi*ht*hr/lambda;    % [m] (alternative definition)
fprintf('Critical distance d_c = %.2f m | Alternative d_break = %.2f m\n', d_c, d_break);
%% -----------------------------------------------------------
% Theoretical models
% -----------------------------------------------------------

% 1) Free-space (Friis) received power in dBm
Lfs_dB = 20*log10(4*pi*d/lambda);  % FSPL
Pr_fs_dBm = Pt_dBm + Gt_dBi + Gr_dBi - Lfs_dB - ap_dB - a_dB;

% 2) Full two-ray with Fresnel ground reflection (complex Γ(θ))
%    Geometry: direct path d1, reflected path d2; grazing angle θ from ground

%Asymptotoic or full two-ray
if max(distVar) < 20*pi*ht*hr/lambda/3
    fprintf('\nDistances are small compared to critical distance; full two-ray model is recommended.\n');
elseif min(distVar) < 20*pi*ht*hr/lambda/3
    fprintf('\nDistances span both near and far regions; full two-ray model is recommended.\n');
end

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

% Create finer distance array for smooth theoretical curve (2x points)
d_fine = linspace(min(d), max(d), 2*length(d));
d1_fine = sqrt(d_fine.^2 + (ht - hr).^2);
d2_fine = sqrt(d_fine.^2 + (ht + hr).^2);
theta_graz_fine = atan2(ht+hr, d_fine);
sin_t_fine = sin(theta_graz_fine);
cos_t_fine = cos(theta_graz_fine);

Xh_fine = sqrt(eps_g - cos_t_fine.^2);
Xv_fine = Xh_fine ./ eps_g;
if upper(pol)=='H'
    X_fine = Xh_fine;
else
    X_fine = Xv_fine;
end
Gamma_fine = (sin_t_fine - X_fine) ./ (sin_t_fine + X_fine);

sumField_fine = exp(-1i*k.*d1_fine)./d1_fine + Gamma_fine .* exp(-1i*k.*d2_fine)./d2_fine;
Pr_2ray_full_W_fine = (10^(Pt_dBm/10)/1000) * Gt * Gr * fac .* (abs(sumField_fine).^2);
Pr_2ray_full_dBm_fine = 10*log10(Pr_2ray_full_W_fine*1000);

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
Pr_2ray_full_dBm_fine_align = Pr_2ray_full_dBm_fine + off_2rf;
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
% Plots: measurement vs two-ray models
% -----------------------------------------------------------
figure(1);
clf;
PL_meas = Pt_dBm + Gt_dBi + Gr_dBi - Rx_meas_dBm;
PL_2ray_full_fine = Pt_dBm + Gt_dBi + Gr_dBi - Pr_2ray_full_dBm_fine_align;
plot(d, PL_meas, 'k-o', 'LineWidth', 3, 'DisplayName','Measured'); hold on;
plot(d_fine, PL_2ray_full_fine, '-', 'LineWidth', 3, 'DisplayName','Two-ray model');
xlabel('Distance (m)', 'FontSize', 32);
ylabel('Path Loss (dB)', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location','best', 'FontSize', 32);
hold off;
saveas(gcf, 'Path_Loss_Two_Ray.png');

% Free-space only comparison
figure(2);
clf;
PL_meas = Pt_dBm + Gt_dBi + Gr_dBi - Rx_meas_dBm;
PL_fs = Pt_dBm + Gt_dBi + Gr_dBi - Pr_fs_dBm_aligned;
plot(d, PL_meas, 'k-o', 'LineWidth', 3, 'DisplayName','Measured'); hold on;
plot(d, PL_fs, '-', 'LineWidth', 3, 'DisplayName',sprintf('Friis (aligned, +%.1f dB)',off_fs));
xlabel('Distance (m)', 'FontSize', 32);
ylabel('Path Loss (dB)', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location','best', 'FontSize', 32);
hold off;
saveas(gcf, 'Free_Space_Path_Loss.png');

%% -----------------------------------------------------------
% Polarization Mismatch Plot
% -----------------------------------------------------------
% Load polarization measurements at 30, 60, and 90 degrees
opts_pol = detectImportOptions('measurements.csv');
opts_pol = setvartype(opts_pol, 'double');
T_pol = readtable('measurements.csv', opts_pol);

% Extract polarization data (columns 3, 4, 5 are 30°, 60°, 90°)
pol_30 = T_pol{:, 3};
pol_60 = T_pol{:, 4};
pol_90 = T_pol{:, 5};
dist_pol = T_pol{:, 1};

% Filter valid measurements for each polarization
valid_30 = ~isnan(pol_30) & ~isnan(dist_pol);
valid_60 = ~isnan(pol_60) & ~isnan(dist_pol);
valid_90 = ~isnan(pol_90) & ~isnan(dist_pol);

figure(3);
clf;

% Define colors for each polarization angle
color_30 = [0.85, 0.33, 0.10];  % Orange-red
color_60 = [0.93, 0.69, 0.13];  % Gold
color_90 = [0.49, 0.18, 0.56];  % Purple

if any(valid_30)
    % Calculate polarization mismatch loss (difference from 0° reference)
    ref_power_30 = interp1(d, Rx_meas_dBm, dist_pol(valid_30), 'linear', 'extrap');
    mismatch_30 = ref_power_30 - pol_30(valid_30);
    
    % Theoretical polarization mismatch loss: L_pol = cos²(θ) in dB = 20*log10(cos(θ))
    theoretical_30 = -20*log10(cosd(30));
    
    plot(dist_pol(valid_30), mismatch_30, 'o-', 'Color', color_30, 'LineWidth', 3, 'DisplayName', '30° Measured');
    hold on;
    yline(theoretical_30, '--', 'Color', color_30, 'DisplayName', sprintf('30° Theoretical (%.2f dB)', theoretical_30), 'LineWidth', 3);
end

if any(valid_60)
    ref_power_60 = interp1(d, Rx_meas_dBm, dist_pol(valid_60), 'linear', 'extrap');
    mismatch_60 = ref_power_60 - pol_60(valid_60);
    theoretical_60 = -20*log10(cosd(60));
    
    plot(dist_pol(valid_60), mismatch_60, 's-', 'Color', color_60, 'LineWidth', 3, 'DisplayName', '60° Measured');
    yline(theoretical_60, '--', 'Color', color_60, 'DisplayName', sprintf('60° Theoretical (%.2f dB)', theoretical_60), 'LineWidth', 3);
end

if any(valid_90)
    ref_power_90 = interp1(d, Rx_meas_dBm, dist_pol(valid_90), 'linear', 'extrap');
    mismatch_90 = ref_power_90 - pol_90(valid_90);
    
    plot(dist_pol(valid_90), mismatch_90, '^-', 'Color', color_90, 'LineWidth', 3, 'DisplayName', '90° Measured');
end

xlabel('Distance (m)', 'FontSize', 32);
ylabel('Polarization Mismatch Loss (dB)', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location','best', 'FontSize', 32);
hold off;
saveas(gcf, 'Polarization_Mismatch_Loss.png');

%% -----------------------------------------------------------
% Antenna Angle Influence Plot
% -----------------------------------------------------------
% Load turning measurements at different angles
opts_turn = detectImportOptions('turning.csv');
opts_turn = setvartype(opts_turn, 'double');
T_turn = readtable('turning.csv', opts_turn);

% Extract angle and power data for both distances
deg_14m = T_turn{:, 1};
rx_14m = T_turn{:, 2};
deg_29m = T_turn{:, 4};
rx_29m = T_turn{:, 5};

% Filter valid measurements
valid_14m = ~isnan(deg_14m) & ~isnan(rx_14m);
valid_29m = ~isnan(deg_29m) & ~isnan(rx_29m);

% Load turntable radiation pattern measurements (with comma as decimal separator)
M_turntable = readmatrix('../turntable/turntable_test.csv', 'DecimalSeparator', ',');

% Extract turntable data
deg_pattern = M_turntable(:, 1);
rx_pattern = M_turntable(:, 2);

% Filter valid turntable measurements and limit to 0-30 degrees range
valid_pattern = ~isnan(deg_pattern) & ~isnan(rx_pattern) & deg_pattern >= 0 & deg_pattern <= 30;
deg_pattern_sub = deg_pattern(valid_pattern);
rx_pattern_sub = rx_pattern(valid_pattern);

% Sort by angle in ascending order for proper plotting
[deg_pattern_sub, sort_idx] = sort(deg_pattern_sub);
rx_pattern_sub = rx_pattern_sub(sort_idx);

% Normalize turntable pattern to 0 dB at 0 degrees (boresight)
idx_zero = find(deg_pattern_sub == 0, 1);
rx_pattern_norm = rx_pattern_sub - rx_pattern_sub(idx_zero);

figure(4);
clf;

% Normalize measurements to 0 dB at 0 degrees for comparison
if any(valid_14m)
    rx_14m_norm = rx_14m(valid_14m) - rx_14m(find(valid_14m, 1, 'first'));
    plot(deg_14m(valid_14m), rx_14m_norm, 'o-', 'LineWidth', 3, 'DisplayName', '14m Measured');
    hold on;
end

if any(valid_29m)
    rx_29m_norm = rx_29m(valid_29m) - rx_29m(find(valid_29m, 1, 'first'));
    plot(deg_29m(valid_29m), rx_29m_norm, 's-', 'LineWidth', 3, 'DisplayName', '29m Measured');
end

% Plot measured radiation pattern from turntable
plot(deg_pattern_sub, rx_pattern_norm, '-', 'LineWidth', 3, 'DisplayName', 'Radiation Pattern (4.5m)');

xlabel('Angle (degrees)', 'FontSize', 32);
ylabel('Relative Power (dB)', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location','best', 'FontSize', 32);
hold off;
saveas(gcf, 'Antenna_Angle_Influence.png');
