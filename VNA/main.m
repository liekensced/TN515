clear all; clc

%% --- S11 (LightGrey) ---
M1 = readmatrix('S11LightGrey.csv', 'Delimiter',';', 'CommentStyle','!', 'DecimalSeparator',',');
f1 = M1(:,1)/1e6;
Gamma1_complex = M1(:,2) + 1i*M1(:,3);  % Complex reflection coefficient
S11_dB = 20*log10(abs(Gamma1_complex));
Gamma1 = abs(Gamma1_complex);  % Magnitude only for VSWR
VSWR1  = (1+Gamma1)./(1-Gamma1);

%% --- S22 (Brick) ---
M2 = readmatrix('S22Brick.csv', 'Delimiter',';', 'CommentStyle','!', 'DecimalSeparator',',');
f2 = M2(:,1)/1e6;
Gamma2_complex = M2(:,2) + 1i*M2(:,3);  % Complex reflection coefficient
S22_dB = 20*log10(abs(Gamma2_complex));
Gamma2 = abs(Gamma2_complex);  % Magnitude only for VSWR
VSWR2  = (1+Gamma2)./(1-Gamma2);

%% --- Plot both one above the other ---
figure('Position',[100 100 900 600]);
set(gca, 'FontSize', 12);

% Top: S11
subplot(2,1,1)
yyaxis left
plot(f1, S11_dB, '-o','LineWidth',1.5,'MarkerFaceColor','b')
ylabel('S11 (dB)'); ylim([-60 0])

yyaxis right
plot(f1, VSWR1, '-s','Color','#D95319','LineWidth',1.5,'MarkerFaceColor','#D95319')
ylabel('VSWR'); ylim([1 10])

xline(3485, 'k-', 'LineWidth', 1.5);
xline((3485+6120)/2, 'k--', 'LineWidth', 1.5);
xline(6120, 'k-', 'LineWidth', 1.5);


grid on; title('S11 - LightGrey Antenna')
xlabel('Frequency (MHz)'); xlim([min(f1) max(f1)])

% Bottom: S22
subplot(2,1,2)
yyaxis left
plot(f2, S22_dB, '-o','LineWidth',1.5,'MarkerFaceColor','b')
ylabel('S22 (dB)'); ylim([-60 0])

yyaxis right
plot(f2, VSWR2, '-s','Color','#D95319','LineWidth',1.5,'MarkerFaceColor','#D95319')
ylabel('VSWR'); ylim([1 10])

xline(3485, 'k-', 'LineWidth', 1.5);
xline((3485+6120)/2, 'k--', 'LineWidth', 1.5);
xline(6120, 'k-', 'LineWidth', 1.5);

grid on; title('S22 - Brick Antenna')
xlabel('Frequency (MHz)'); xlim([min(f2) max(f2)])

%% --- Reflection Coefficient Plot ---
figure('Position',[1050 100 900 600]);
set(gca, 'FontSize', 12);
plot(f1, Gamma1, '-o','LineWidth',1.5,'MarkerFaceColor','b', 'DisplayName', 'LightGrey (S11)')
hold on
plot(f2, Gamma2, '-s','LineWidth',1.5,'MarkerFaceColor','r', 'DisplayName', 'Brick (S22)')
ylabel('|\Gamma|'); xlabel('Frequency (MHz)');
grid on; title('Reflection Coefficient Magnitude')
legend('Location','best')
xlim([min(min(f1),min(f2)) max(max(f1),max(f2))])
ylim([0 1])

xline(3485, 'k-', 'LineWidth', 1.5);
xline((3485+6120)/2, 'k--', 'LineWidth', 1.5);
xline(6120, 'k-', 'LineWidth', 1.5);

saveas(gcf, 'S11_S22_Comparison.png');

%% --- Calculate metrics in working domain ---
f_min = 3485;  % MHz
f_max = 6120;  % MHz
f_center = (f_min + f_max) / 2;  % 4802.5 MHz

% Find indices in working domain for both antennas
idx1_working = f1 >= f_min & f1 <= f_max;
idx2_working = f2 >= f_min & f2 <= f_max;

% LightGrey Antenna (S11)
fprintf('\n========== LightGrey Antenna (S11) ==========\n');
fprintf('Working Domain: %.0f - %.0f MHz\n\n', f_min, f_max);

% Find closest frequency to center
[~, idx1_center] = min(abs(f1 - f_center));
fprintf('At Center Frequency (%.1f MHz):\n', f1(idx1_center));
fprintf('  Return Loss (dB):    %.2f\n', S11_dB(idx1_center));
fprintf('  VSWR:                %.4f\n', VSWR1(idx1_center));
fprintf('  Reflection Coeff:    %.4f\n', Gamma1(idx1_center));

% In working domain
fprintf('\nIn Working Domain:\n');
fprintf('  Max Return Loss (dB):    %.2f\n', max(S11_dB(idx1_working)));
fprintf('  Max VSWR:                %.4f\n', max(VSWR1(idx1_working)));
fprintf('  Max Reflection Coeff:    %.4f\n', max(Gamma1(idx1_working)));
fprintf('  Avg Return Loss (dB):    %.2f\n', mean(S11_dB(idx1_working)));
fprintf('  Avg VSWR:                %.4f\n', mean(VSWR1(idx1_working)));
fprintf('  Avg Reflection Coeff:    %.4f\n', mean(Gamma1(idx1_working)));

% Brick Antenna (S22)
fprintf('\n========== Brick Antenna (S22) ==========\n');
fprintf('Working Domain: %.0f - %.0f MHz\n\n', f_min, f_max);

% Find closest frequency to center
[~, idx2_center] = min(abs(f2 - f_center));
fprintf('At Center Frequency (%.1f MHz):\n', f2(idx2_center));
fprintf('  Return Loss (dB):    %.2f\n', S22_dB(idx2_center));
fprintf('  VSWR:                %.4f\n', VSWR2(idx2_center));
fprintf('  Reflection Coeff:    %.4f\n', Gamma2(idx2_center));

% In working domain
fprintf('\nIn Working Domain:\n');
fprintf('  Max Return Loss (dB):    %.2f\n', max(S22_dB(idx2_working)));
fprintf('  Max VSWR:                %.4f\n', max(VSWR2(idx2_working)));
fprintf('  Max Reflection Coeff:    %.4f\n', max(Gamma2(idx2_working)));
fprintf('  Avg Return Loss (dB):    %.2f\n', mean(S22_dB(idx2_working)));
fprintf('  Avg VSWR:                %.4f\n', mean(VSWR2(idx2_working)));
fprintf('  Avg Reflection Coeff:    %.4f\n\n', mean(Gamma2(idx2_working)));

% Save the figure
saveas(gcf, 'Reflection.png');

%% --- Impedance Matching Analysis at Center Frequency ---
fprintf('========== Impedance Matching Analysis (f = %.1f MHz) ==========\n', f_center);
Z0 = 50;  % Characteristic impedance in Ohms

% LightGrey Antenna
Gamma1_center = Gamma1_complex(idx1_center);  % Use complex value
Z_load1 = Z0 * (1 + Gamma1_center) / (1 - Gamma1_center);
Z_load1_real = real(Z_load1);
Z_load1_imag = imag(Z_load1);
Z_qw1 = sqrt(Z0 * abs(Z_load1));

fprintf('\nLightGrey Antenna (S11):\n');
fprintf('  Antenna Impedance:   %.2f + j%.2f Ω\n', Z_load1_real, Z_load1_imag);
fprintf('  λ/4 Transformer Z:   %.2f Ω\n', Z_qw1);

% Brick Antenna
Gamma2_center = Gamma2_complex(idx2_center);  % Use complex value
Z_load2 = Z0 * (1 + Gamma2_center) / (1 - Gamma2_center);
Z_load2_real = real(Z_load2);
Z_load2_imag = imag(Z_load2);
Z_qw2 = sqrt(Z0 * abs(Z_load2));

fprintf('\nBrick Antenna (S22):\n');
fprintf('  Antenna Impedance:   %.2f + j%.2f Ω\n', Z_load2_real, Z_load2_imag);
fprintf('  λ/4 Transformer Z:   %.2f Ω\n\n', Z_qw2);




%sgtitle('S11 vs S22 Comparison   (dual axis: dB left, VSWR right)')