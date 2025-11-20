clear; clc

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