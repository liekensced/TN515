clear; clc; close all;

%% Horn Antenna Simulation at 4.8 GHz
% This script simulates a pyramidal horn antenna using MATLAB Antenna Toolbox
% Dimensions (converted to meters):
% - Aperture width (FlareWidth): 10.4 cm = 0.104 m
% - Aperture height (FlareHeight): 13.8 cm = 0.138 m
% - Inner rib: 18.5 cm = 0.185 m -> Flare length = 0.174 m
% - Throat/waveguide width (Width): 2.2 cm = 0.022 m
% - Throat/waveguide height (Height
%   ): 4.75 cm = 0.0475 m
% Frequency: 4.8 GHz

% Create default horn antenna
h = horn;
h.FlareWidth = 0.138; h.FlareHeight = 0.104; h.FlareLength = 0.174;
h.Width = 0.0475; h.Height = 0.022; h.Length = 0.05;
freq = 4.8e9;  % 4.8 GHz

% Display antenna info
figure;
show(h);
title('Pyramidal Horn Antenna Geometry');

% Plot 3D radiation pattern
figure;
pattern(h, freq);
title('3D Radiation Pattern at 4.8 GHz');
[pat, az_3d, el_3d] = pattern(h, freq);
Gpeak = max(pat(:));
fprintf('Simulated Peak Directivity: %.2f dBi\n', Gpeak);

% Extract azimuth (H-plane) pattern data
az = -180:180;
[pat_az] = patternAzimuth(h, freq, 0, 'Azimuth', az);

% Extract elevation (E-plane) pattern data
el = -180:180;
[pat_el] = patternElevation(h, freq, 90, 'Elevation', el);


% Far-field approximation distance
D = max(h.FlareWidth, h.FlareHeight); % Max dimension
lambda = physconst('LightSpeed') / freq;
R_far = 2 * D^2 / lambda;
fprintf('Recommended far-field distance: %.2f m\n', R_far);

%% Test Data Processing
% 1) Azimuth pattern (turntable_test.csv)
M1 = readmatrix('turntable_test.csv', 'DecimalSeparator',',');
az_deg = M1(:,1);
az_dBm = M1(:,2);
az_norm = az_dBm - max(az_dBm);

% 2) Elevation pattern (turntable_test_elevation.csv)
M2 = readmatrix('turntable_test_elevation.csv', 'DecimalSeparator',',');
el_deg = M2(:,1);
el_raw = M2(:,2);
el_dBm = -el_raw; % fix the sign!
el_norm = el_dBm - max(el_dBm);

%% Parameters for gain calculation
f = 4.8e9; % 4.8 GHz
c = 3e8;
lambda = c/f; % wavelength = 0.0625 m
d = 4.5; % distance = 4.5 m
Pt_dBm = 10; % transmitted power = 10 dBm
Pt = 10^((Pt_dBm-30)/10); % Pt in Watts
% Free-space path loss (ideal, L=1)
FSPL = (lambda / (4*pi*d))^2; % linear
FSPL_dB = 10*log10(FSPL); % dB
% Theoretical received power in boresight (max gain direction) if isotropic antennas
Pr_iso_dBm = Pt_dBm + 20*log10(lambda/(4*pi*d)); % = Pt_dBm - FSPL_dB

%% Find realized gain in boresight (maximum measured direction)
% Measured max received power (in dBm)
Pr_meas_max_dBm_az = max(az_dBm);
Pr_meas_max_dBm_el = max(el_dBm);
% Since both antennas are identical → Gt = Gr = G
% Friis: Pr = Pt * Gt * Gr * (lambda/(4*pi*d))^2 → Pr/Pt = G^2 * FSPL
G_linear_az = sqrt( 10.^((Pr_meas_max_dBm_az - Pt_dBm)/10) / FSPL );
G_linear_el = sqrt( 10.^((Pr_meas_max_dBm_el - Pt_dBm)/10) / FSPL );
G_dBi_az = 10*log10(G_linear_az);
G_dBi_el = 10*log10(G_linear_el);
fprintf('--- Antenna Gain Calculation (4.8 GHz, d = 4.5 m) ---\n');
fprintf('Free-space path loss (FSPL) = %.2f dB\n', -FSPL_dB);
fprintf('Theoretical Pr (isotropic) = %.2f dBm\n', Pr_iso_dBm);
fprintf('Measured max Pr (azimuth cut) = %.2f dBm\n', Pr_meas_max_dBm_az);
fprintf('Measured max Pr (elevation cut) = %.2f dBm\n', Pr_meas_max_dBm_el);
fprintf('Realized Gain (azimuth boresight) = %.2f dBi\n', G_dBi_az);
fprintf('Realized Gain (elevation boresight) = %.2f dBi\n', G_dBi_el);
fprintf('Average boresight gain = %.2f dBi\n', mean([G_dBi_az G_dBi_el]));

%% Convert normalized patterns to absolute gain patterns (dBi)
az_gain_dBi = az_norm + max(G_dBi_az, G_dBi_el); % use the higher of the two cuts
el_gain_dBi = el_norm + max(G_dBi_az, G_dBi_el);

%% Plot both simulated and test patterns
figure('Position',[100 100 1000 450])
% Azimuth
subplot(1,2,1)
polarplot(deg2rad(az_deg), az_gain_dBi, '-o','LineWidth',2,'MarkerFaceColor','b','DisplayName','Test')
hold on
polarplot(deg2rad(az), pat_az, 'r--','LineWidth',2,'DisplayName','Simulation')
title('Azimuth Pattern (H-plane)')
rlim([-40 max([0, max(az_gain_dBi), max(pat_az)])])
rticks(-40:10:max([0, max(az_gain_dBi), max(pat_az)]))
thetaticks(0:45:315)
ax = gca; ax.ThetaZeroLocation = 'right';
legend('Location','best')

% Elevation
subplot(1,2,2)
polarplot(deg2rad(el_deg-90), el_gain_dBi, '-o','LineWidth',2,'MarkerFaceColor','r','DisplayName','Test')
hold on
polarplot(deg2rad(el), pat_el, 'b--','LineWidth',2,'DisplayName','Simulation')
title('Elevation Pattern (E-plane) - VERTICAL')
rlim([-40 max([0, max(el_gain_dBi), max(pat_el)])])
rticks(-40:10:max([0, max(el_gain_dBi), max(pat_el)]))
thetaticks(-90:30:90)
ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise';
legend('Location','best')

sgtitle('4.8 GHz Radiation Patterns: Simulation vs Test (4.5 m, 10 dBm TX)')

%% Plot normalized patterns (relative to each cut's own maximum)
figure('Position',[100 100 1000 450])
% Azimuth - normalized
subplot(1,2,1)
az_gain_norm = az_norm; % already normalized to 0 dB max
pat_az_norm = pat_az - max(pat_az); % normalize simulation
polarplot(deg2rad(az_deg), az_gain_norm, '-o','LineWidth',2,'MarkerFaceColor','b','DisplayName','Test')
hold on
polarplot(deg2rad(az), pat_az_norm, 'r--','LineWidth',2,'DisplayName','Simulation')
title('Azimuth Pattern (H-plane)')
rlim([-40 0])
rticks(-40:10:0)
thetaticks(0:45:315)
ax = gca; ax.ThetaZeroLocation = 'right';
legend('Location','best')

% Elevation - normalized
subplot(1,2,2)
el_gain_norm = el_norm; % already normalized to 0 dB max
pat_el_norm = pat_el - max(pat_el); % normalize simulation
polarplot(deg2rad(el_deg-90), el_gain_norm, '-o','LineWidth',2,'MarkerFaceColor','r','DisplayName','Test')
hold on
polarplot(deg2rad(el), pat_el_norm, 'b--','LineWidth',2,'DisplayName','Simulation')
title('Elevation Pattern (E-plane) - VERTICAL')
rlim([-40 0])
rticks(-40:10:0)
thetaticks(-90:30:90)
ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise';
legend('Location','best')

sgtitle('4.8 GHz Radiation Patterns: Normalized (0 dB at peak)')

%% Key Metrics Calculation
fprintf('\n=== Key Metrics Report ===\n');

%% Azimuth (H-plane) Test
[beamwidth_az_test, sll_az_test, fb_az_test, null_depths_az_test, null_locs_az_test, peak_az_test] = compute_metrics(az_deg, az_gain_dBi);

%% Azimuth (H-plane) Simulation
[beamwidth_az_sim, sll_az_sim, fb_az_sim, null_depths_az_sim, null_locs_az_sim, peak_az_sim] = compute_metrics(az, pat_az);

%% Elevation (E-plane) Test
[beamwidth_el_test, sll_el_test, fb_el_test, null_depths_el_test, null_locs_el_test, peak_el_test] = compute_metrics(el_deg, el_gain_dBi);

%% Elevation (E-plane) Simulation
[beamwidth_el_sim, sll_el_sim, fb_el_sim, null_depths_el_sim, null_locs_el_sim, peak_el_sim] = compute_metrics(el, pat_el);

%% Pointing Error
pointing_error_az = min(abs(peak_az_test - peak_az_sim), 360 - abs(peak_az_test - peak_az_sim));
pointing_error_el = min(abs(peak_el_test - peak_el_sim), 360 - abs(peak_el_test - peak_el_sim));

%% Report
fprintf('Azimuth (H-plane):\n');
fprintf('  Test: 3dB BW=%.1f°, SLL=%.1f dB, F/B=%.1f dB, Peak at %.1f°\n', beamwidth_az_test, sll_az_test, fb_az_test, peak_az_test);
fprintf('  Sim:  3dB BW=%.1f°, SLL=%.1f dB, F/B=%.1f dB, Peak at %.1f°\n', beamwidth_az_sim, sll_az_sim, fb_az_sim, peak_az_sim);
fprintf('  Pointing Error: %.1f°\n', pointing_error_az);
fprintf('  Null Depths/Locs (Test): ');
for i = 1:min(3, length(null_depths_az_test))
    fprintf('%.1f dB at %.1f°, ', null_depths_az_test(i), null_locs_az_test(i));
end
fprintf('\n');
fprintf('  Null Depths/Locs (Sim): ');
for i = 1:min(3, length(null_depths_az_sim))
    fprintf('%.1f dB at %.1f°, ', null_depths_az_sim(i), null_locs_az_sim(i));
end
fprintf('\n');

fprintf('Elevation (E-plane):\n');
fprintf('  Test: 3dB BW=%.1f°, SLL=%.1f dB, F/B=%.1f dB, Peak at %.1f°\n', beamwidth_el_test, sll_el_test, fb_el_test, peak_el_test);
fprintf('  Sim:  3dB BW=%.1f°, SLL=%.1f dB, F/B=%.1f dB, Peak at %.1f°\n', beamwidth_el_sim, sll_el_sim, fb_el_sim, peak_el_sim);
fprintf('  Pointing Error: %.1f°\n', pointing_error_el);
fprintf('  Null Depths/Locs (Test): ');
for i = 1:min(3, length(null_depths_el_test))
    fprintf('%.1f dB at %.1f°, ', null_depths_el_test(i), null_locs_el_test(i));
end
fprintf('\n');
fprintf('  Null Depths/Locs (Sim): ');
for i = 1:min(3, length(null_depths_el_sim))
    fprintf('%.1f dB at %.1f°, ', null_depths_el_sim(i), null_locs_el_sim(i));
end
fprintf('\n');

%% Helper function to compute metrics

function [beamwidth, sll, fb_ratio_db, null_depths, null_locs, peak_ang] = compute_metrics(angles_deg, gains_dB)
% COMPUTE_METRICS Robust antenna pattern metrics on circular angle grids (degrees).
% Inputs:
%   angles_deg : angle vector (deg), can be 0..360, -180..180, etc. Order/duplicates OK.
%   gains_dB   : gain/pattern (dB) at each angle.
%
% Outputs:
%   beamwidth     : 3 dB beamwidth (deg), NaN if cannot be determined
%   sll           : first sidelobe level RELATIVE to peak (<= 0 dB)
%   fb_ratio_db   : front-to-back ratio (hemispheric), peak - max(back) (>= 0 dB)
%   null_depths   : null depths RELATIVE to peak (<= 0 dB)
%   null_locs     : null locations (deg), normalized to [-180, 180)
%   peak_ang      : peak direction (deg), normalized to [-180, 180)
%
% Notes:
% - Angles are normalized to [-180, 180) internally, duplicates removed, then sorted.
% - 3 dB beamwidth is found via circular search and linear interpolation of threshold crossings.
% - SLL is taken between the first two nulls to the right of the main beam (circularly).
% - F/B is hemispheric: max gain in [peak+90°, peak+270°], then peak - back_max.

% ----------- Input prep: vectorize, clean, normalize, dedupe, sort -----------
a = angles_deg(:);
g = gains_dB(:);

% Remove non-finite samples
finiteMask = isfinite(a) & isfinite(g);
a = a(finiteMask);
g = g(finiteMask);

if numel(a) < 3
    beamwidth   = NaN;
    sll         = NaN;
    fb_ratio_db = NaN;
    null_depths = NaN(0,1);
    null_locs   = NaN(0,1);
    peak_ang    = NaN;
    return
end

% Normalize angles to [-180, 180) and remove duplicates after normalization
a = mod(a + 180, 360) - 180;
[a, uniq] = unique(a, 'stable');
g = g(uniq);

% Sort by angle
[a, sidx] = sort(a);
g = g(sidx);
N = numel(a);

% Guard: ensure we still have enough samples
if N < 3
    beamwidth   = NaN;
    sll         = NaN;
    fb_ratio_db = NaN;
    null_depths = NaN(0,1);
    null_locs   = NaN(0,1);
    peak_ang    = NaN;
    return
end

% ----------- Peak detection -----------
[peak_gain, idx_peak] = max(g);
peak_ang = a(idx_peak);  % will be normalized later anyway

% ----------- 3 dB beamwidth (circular with interpolation) -----------
half_power = peak_gain - 3;

% Create triple-copy arrays for safe wrap-around indexing
A3 = [a - 360; a; a + 360];
G3 = [g; g; g];
ip = idx_peak + N;  % peak index in central copy
M  = numel(G3);

% If pattern never drops below half-power (near-omni), beamwidth undefined
if ~any(G3 < half_power)
    beamwidth = NaN;
else
    % ---- Left crossing (moving left from peak until we drop below threshold) ----
    j = ip;
    while j > 1 && G3(j) >= half_power
        j = j - 1;
    end
    left_ok = (j >= 1) && (j < M);
    if left_ok
        % Interpolate between j (below) and j+1 (at/above) to get the crossing
        if G3(j) == G3(j+1)
            left3dB = A3(j); % flat segment on threshold (rare); pick left point
        else
            left3dB = interp1(G3(j:j+1), A3(j:j+1), half_power, 'linear', 'extrap');
        end
    else
        left3dB = NaN;
    end
    
    % ---- Right crossing (moving right from peak until we drop below threshold) ----
    k = ip;
    while k < M && G3(k) >= half_power
        k = k + 1;
    end
    right_ok = (k > 1) && (k <= M);
    if right_ok
        % Interpolate between k-1 (at/above) and k (below)
        if G3(k-1) == G3(k)
            right3dB = A3(k); % flat segment on threshold (rare)
        else
            right3dB = interp1(G3(k-1:k), A3(k-1:k), half_power, 'linear', 'extrap');
        end
    else
        right3dB = NaN;
    end
    
    if isnan(left3dB) || isnan(right3dB)
        beamwidth = NaN;
    else
        beamwidth = right3dB - left3dB;
        % sanity: if interpolation/ordering hiccups produce non-sensical width
        if beamwidth <= 0 || beamwidth >= 360
            beamwidth = NaN;
        end
    end
end

% ----------- Nulls (local minima) and relative depths -----------
% Find local minima on triple array
[negpks_all, locs_all] = findpeaks(-G3);  % minima of G3
null_vals_all = -negpks_all;
null_angs_all = A3(locs_all);

% Keep only nulls within the central period to avoid duplicates
in_central = locs_all > N & locs_all <= 2*N;
null_vals  = null_vals_all(in_central);
null_angs  = null_angs_all(in_central);

% Sort by angle
[null_angs, nsi] = sort(null_angs);
null_vals = null_vals(nsi);

% Relative null depths (<= 0 dB), locations normalized
null_depths = null_vals - peak_gain;
null_locs   = mod(null_angs + 180, 360) - 180;

% ----------- First sidelobe level (relative) -----------
% Use nulls to the right of the peak (circularly) from the triple array to be safe
right_mask_all = null_angs_all > peak_ang;              % in A3's angle space
if ~any(right_mask_all)
    sll = NaN;
else
    % Consider nulls to the right within one full turn
    ra = null_angs_all(right_mask_all);
    rv = null_vals_all(right_mask_all);
    % Keep those within (peak, peak+360]
    sel = ra <= (peak_ang + 360);
    ra = ra(sel); rv = rv(sel);
    [ra, rsi] = sort(ra); rv = rv(rsi);
    
    if numel(ra) < 2
        sll = NaN;
    else
        first_null  = ra(1);
        second_null = ra(2);
        
        % Interval strictly between first two right-hand nulls
        mask = (A3 > first_null) & (A3 < second_null);
        if any(mask)
            sidelobe_max = max(G3(mask));
            sll = sidelobe_max - peak_gain;  % <= 0 dB
        else
            sll = NaN;
        end
    end
end

% ----------- Front-to-Back ratio (hemispheric) -----------
back_lo = peak_ang + 90;
back_hi = peak_ang + 270;
back_mask = (A3 >= back_lo) & (A3 <= back_hi);
if any(back_mask)
    back_max = max(G3(back_mask));
    fb_ratio_db = peak_gain - back_max;   % >= 0 dB
else
    fb_ratio_db = NaN;
end

% ----------- Normalize output angles to [-180, 180) -----------
peak_ang = mod(peak_ang + 180, 360) - 180;

end