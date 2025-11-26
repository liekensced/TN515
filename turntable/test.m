clear; clc
import simulation
%% 1) Azimuth pattern (turntable_test.csv)
M1 = readmatrix('turntable_test.csv', 'DecimalSeparator',',');
az_deg = M1(:,1);
az_dBm = M1(:,2);
az_norm = az_dBm - max(az_dBm);

%% 2) Elevation pattern (turntable_test_elevation.csv)
M2 = readmatrix('turntable_test_elevation.csv', 'DecimalSeparator',',');
el_deg = M2(:,1);
el_raw = M2(:,2);
el_dBm = -el_raw;                    % fix the sign!
el_norm = el_dBm - max(el_dBm);

%% Plot both correctly (no yyaxis, no errors)
figure('Position',[100 100 1000 450])

% Azimuth
subplot(1,2,1)
polarplot(deg2rad(az_deg), az_norm, '-o','LineWidth',2,'MarkerFaceColor','b')
title('Azimuth Pattern (H-plane)')
rlim([-40 0])
rticks(-40:10:0)
thetaticks(0:45:315)
ax = gca; ax.ThetaZeroLocation = 'right';

% Elevation
subplot(1,2,2)
polarplot(deg2rad(el_deg-90), el_norm, '-o','LineWidth',2,'MarkerFaceColor','r')
title('Elevation Pattern (E-plane) - VERTICAL')
rlim([-40 0])
rticks(-40:10:0)
thetaticks(-90:30:90)
ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise';

sgtitle('4.8 GHz Radiation Patterns  (4.5 m, 10 dBm TX)')

%% Parameters for gain calculation
f  = 4.8e9;          % 4.8 GHz
c  = 3e8;
lambda = c/f;         % wavelength = 0.0625 m
d  = 4.5;             % distance = 4.5 m
Pt_dBm = 10;          % transmitted power = 10 dBm
Pt = 10^((Pt_dBm-30)/10);  % Pt in Watts

% Free-space path loss (ideal, L=1)
FSPL = (lambda / (4*pi*d))^2;           % linear
FSPL_dB = 10*log10(FSPL);               % dB

% Theoretical received power in boresight (max gain direction) if isotropic antennas
Pr_iso_dBm = Pt_dBm + 20*log10(lambda/(4*pi*d));   % = Pt_dBm - FSPL_dB

%% Find realized gain in boresight (maximum measured direction)
% Measured max received power (in dBm)
Pr_meas_max_dBm_az = max(az_dBm);
Pr_meas_max_dBm_el = max(el_dBm);

% Since both antennas are identical → Gt = Gr = G
% Friis: Pr = Pt * Gt * Gr * (lambda/(4*pi*d))^2  →  Pr/Pt = G^2 * FSPL
G_linear_az = sqrt( 10.^((Pr_meas_max_dBm_az - Pt_dBm)/10) / FSPL );
G_linear_el = sqrt( 10.^((Pr_meas_max_dBm_el - Pt_dBm)/10) / FSPL );

G_dBi_az = 10*log10(G_linear_az);
G_dBi_el = 10*log10(G_linear_el);

fprintf('--- Antenna Gain Calculation (4.8 GHz, d = 4.5 m) ---\n');
fprintf('Free-space path loss (FSPL)      = %.2f dB\n', -FSPL_dB);
fprintf('Theoretical Pr (isotropic)      = %.2f dBm\n', Pr_iso_dBm);
fprintf('Measured max Pr (azimuth cut)    = %.2f dBm\n', Pr_meas_max_dBm_az);
fprintf('Measured max Pr (elevation cut)  = %.2f dBm\n', Pr_meas_max_dBm_el);
fprintf('Realized Gain (azimuth boresight)   = %.2f dBi\n', G_dBi_az);
fprintf('Realized Gain (elevation boresight) = %.2f dBi\n', G_dBi_el);
fprintf('Average boresight gain              = %.2f dBi\n', mean([G_dBi_az G_dBi_el]));

%% Convert normalized patterns to absolute gain patterns (dBi)
az_gain_dBi = az_norm + max(G_dBi_az, G_dBi_el);   % use the higher of the two cuts
el_gain_dBi = el_norm + max(G_dBi_az, G_dBi_el);