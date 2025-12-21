% Plot fading dBm samples produced by GNU Radio
% Assumes 32-bit float little-endian samples (dBm) written raw to file.

% Get the directory where this script is located
scriptDir = fileparts(mfilename('fullpath'));

sampleRateHz = 1.28e6;           % capture sample speed (Hz)
fftSize = 8192;                  % FFT size used in GNU Radio
bytesPerSample = 4;              % float32
frameRate = sampleRateHz / fftSize;  % frames per second (power updates)

% Load LOS file
fileName_LOS = 'LOS_2048_Sam.txt';
fid = fopen(fileName_LOS, 'rb');
if fid < 0
    error('Could not open file: %s', fileName_LOS);
end
samples_LOS = fread(fid, inf, 'float32', 0, 'l');
fclose(fid);
timeVec_LOS = (0:numel(samples_LOS)-1).' / frameRate;

% Load NOLOS file
fileName_NOLOS = 'NOLOS_2048_Sam.txt';
fid = fopen(fileName_NOLOS, 'rb');
if fid < 0
    error('Could not open file: %s', fileName_NOLOS);
end
samples_NOLOS = fread(fid, inf, 'float32', 0, 'l');
fclose(fid);
timeVec_NOLOS = (0:numel(samples_NOLOS)-1).' / frameRate;

% Plot LOS
figure(1);
clf;
plot(timeVec_LOS, samples_LOS, '-', 'LineWidth', 3);
xlabel('Time (s)', 'FontSize', 32);
ylabel('Rx (dBm)', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
saveas(gcf, fullfile(scriptDir, 'LOS_2048.png'));

% Plot NLOS
figure(2);
clf;
plot(timeVec_NOLOS, samples_NOLOS, '-', 'LineWidth', 3);
xlabel('Time (s)', 'FontSize', 32);
ylabel('Rx (dBm)', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
saveas(gcf, fullfile(scriptDir, 'NLOS_2048.png'));

%% -----------------------------------------------------------
% Distribution Analysis
% -----------------------------------------------------------

% Load all NOLOS files
fid = fopen('NOLOS_512_Sam.txt', 'rb');
samples_NOLOS_512 = fread(fid, inf, 'float32', 0, 'l');
fclose(fid);

fid = fopen('NOLOS_8192_Sam.txt', 'rb');
samples_NOLOS_8192 = fread(fid, inf, 'float32', 0, 'l');
fclose(fid);

% Load all LOS files
fid = fopen('LOS_512_Sam.txt', 'rb');
samples_LOS_512 = fread(fid, inf, 'float32', 0, 'l');
fclose(fid);

fid = fopen('LOS_8192_Sam.txt', 'rb');
samples_LOS_8192 = fread(fid, inf, 'float32', 0, 'l');
fclose(fid);

% Convert dBm to linear voltage (amplitude) for distribution fitting
% P(mW) = 10^(P_dBm/10), then V ~ sqrt(P) for amplitude
dBm_to_amp = @(x) sqrt(10.^(x/10));

amp_NOLOS_512 = dBm_to_amp(samples_NOLOS_512);
amp_NOLOS_2048 = dBm_to_amp(samples_NOLOS);
amp_NOLOS_8192 = dBm_to_amp(samples_NOLOS_8192);

amp_LOS_512 = dBm_to_amp(samples_LOS_512);
amp_LOS_2048 = dBm_to_amp(samples_LOS);
amp_LOS_8192 = dBm_to_amp(samples_LOS_8192);

% NOLOS - Rayleigh Distribution Comparison
figure(3);
clf;
hold on;

% Define colors for each FFT size
color_512 = [0.5, 0, 0.50];   % Purple
color_2048 = [1, 0, 0.];  % Red
color_8192 = [0.00, 0, 1];  % Blue

% Calculate RMS (sigma) for each NOLOS dataset using MLE for Rayleigh
% For Rayleigh: B = sqrt(mean(r^2)/2)
sigma_NOLOS_512 = sqrt(mean(amp_NOLOS_512.^2) / 2);
sigma_NOLOS_2048 = sqrt(mean(amp_NOLOS_2048.^2) / 2);
sigma_NOLOS_8192 = sqrt(mean(amp_NOLOS_8192.^2) / 2);

% Rayleigh PDF: r/sigma^2 * exp(-r^2/(2*sigma^2))
rayleigh_pdf = @(r, sigma) (r ./ sigma^2) .* exp(-r.^2 ./ (2*sigma^2));

% Create histogram and overlay theoretical
[counts_512, edges_512] = histcounts(amp_NOLOS_512, 50, 'Normalization', 'pdf');
centers_512 = (edges_512(1:end-1) + edges_512(2:end))/2;
plot(centers_512, counts_512, 'o', 'Color', color_512, 'LineWidth', 3, 'DisplayName', 'NLOS 512');

[counts_2048, edges_2048] = histcounts(amp_NOLOS_2048, 50, 'Normalization', 'pdf');
centers_2048 = (edges_2048(1:end-1) + edges_2048(2:end))/2;
plot(centers_2048, counts_2048, 's', 'Color', color_2048, 'LineWidth', 3, 'DisplayName', 'NLOS 2048');

[counts_8192, edges_8192] = histcounts(amp_NOLOS_8192, 50, 'Normalization', 'pdf');
centers_8192 = (edges_8192(1:end-1) + edges_8192(2:end))/2;
plot(centers_8192, counts_8192, '^', 'Color', color_8192, 'LineWidth', 3, 'DisplayName', 'NLOS 8192');

% Plot theoretical Rayleigh distributions for each dataset
r_ray = linspace(0, max([amp_NOLOS_512; amp_NOLOS_2048; amp_NOLOS_8192]), 200);
plot(r_ray, rayleigh_pdf(r_ray, sigma_NOLOS_512), '--', 'Color', color_512, 'LineWidth', 3, 'DisplayName', 'Rayleigh 512');
plot(r_ray, rayleigh_pdf(r_ray, sigma_NOLOS_2048), '--', 'Color', color_2048, 'LineWidth', 3, 'DisplayName', 'Rayleigh 2048');
plot(r_ray, rayleigh_pdf(r_ray, sigma_NOLOS_8192), '--', 'Color', color_8192, 'LineWidth', 3, 'DisplayName', 'Rayleigh 8192');

xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'NLOS_Rayleigh_Distribution.png'));

% Goodness-of-fit test: Kolmogorov-Smirnov test for Rayleigh
fprintf('\n--- NLOS Rayleigh Distribution K-S Test ---\n');
[h_512, p_512, ks_512] = kstest(amp_NOLOS_512, 'CDF', makedist('Rayleigh', 'B', sigma_NOLOS_512));
fprintf('NLOS 512:  K-S statistic = %.4f, p-value = %.4f %s\n', ks_512, p_512, fit_quality(ks_512));

[h_2048, p_2048, ks_2048] = kstest(amp_NOLOS_2048, 'CDF', makedist('Rayleigh', 'B', sigma_NOLOS_2048));
fprintf('NLOS 2048: K-S statistic = %.4f, p-value = %.4f %s\n', ks_2048, p_2048, fit_quality(ks_2048));

[h_8192, p_8192, ks_8192] = kstest(amp_NOLOS_8192, 'CDF', makedist('Rayleigh', 'B', sigma_NOLOS_8192));
fprintf('NLOS 8192: K-S statistic = %.4f, p-value = %.4f %s\n', ks_8192, p_8192, fit_quality(ks_8192));
fprintf('Note: K-S statistic interpretation: <0.05=excellent, <0.10=good, <0.20=reasonable\n');

% LOS - Rician Distribution Comparison
figure(4);
clf;
hold on;

% Define colors for each FFT size (same as NLOS)
color_512 = [0.5, 0, 0.50];   % Purple
color_2048 = [1, 0, 0.];  % Red
color_8192 = [0.00, 0, 1];  % Blue

% Calculate parameters for each LOS dataset
A_LOS_512 = mean(amp_LOS_512);
sigma_LOS_512 = std(amp_LOS_512);

A_LOS_2048 = mean(amp_LOS_2048);
sigma_LOS_2048 = std(amp_LOS_2048);

A_LOS_8192 = mean(amp_LOS_8192);
sigma_LOS_8192 = std(amp_LOS_8192);

% Calculate Rician K-factor (in dB) for each dataset
% K = A^2 / (2*sigma^2)
K_dB_512 = 10*log10(A_LOS_512^2 / (2*sigma_LOS_512^2));
K_dB_2048 = 10*log10(A_LOS_2048^2 / (2*sigma_LOS_2048^2));
K_dB_8192 = 10*log10(A_LOS_8192^2 / (2*sigma_LOS_8192^2));

% Rician PDF: r/sigma^2 * exp(-(r^2 + A^2)/(2*sigma^2)) * I_0(A*r/sigma^2)
% I_0 is modified Bessel function of first kind, zero order
rician_pdf = @(r, A, sigma) (r ./ sigma^2) .* exp(-(r.^2 + A^2) ./ (2*sigma^2)) .* besseli(0, A*r ./ sigma^2);

% Create histogram and overlay theoretical
[counts_los_512, edges_los_512] = histcounts(amp_LOS_512, 50, 'Normalization', 'pdf');
centers_los_512 = (edges_los_512(1:end-1) + edges_los_512(2:end))/2;
plot(centers_los_512, counts_los_512, 'o', 'Color', color_512, 'LineWidth', 3, 'DisplayName', 'LOS 512');

[counts_los_2048, edges_los_2048] = histcounts(amp_LOS_2048, 50, 'Normalization', 'pdf');
centers_los_2048 = (edges_los_2048(1:end-1) + edges_los_2048(2:end))/2;
plot(centers_los_2048, counts_los_2048, 's', 'Color', color_2048, 'LineWidth', 3, 'DisplayName', 'LOS 2048');

[counts_los_8192, edges_los_8192] = histcounts(amp_LOS_8192, 50, 'Normalization', 'pdf');
centers_los_8192 = (edges_los_8192(1:end-1) + edges_los_8192(2:end))/2;
plot(centers_los_8192, counts_los_8192, '^', 'Color', color_8192, 'LineWidth', 3, 'DisplayName', 'LOS 8192');

% Plot theoretical Rician distributions for each dataset
r_rice = linspace(0, max([amp_LOS_512; amp_LOS_2048; amp_LOS_8192]), 200);
plot(r_rice, rician_pdf(r_rice, A_LOS_512, sigma_LOS_512), '--', 'Color', color_512, 'LineWidth', 3, 'DisplayName', sprintf('Rician 512 (K=%.1f dB)', K_dB_512));
plot(r_rice, rician_pdf(r_rice, A_LOS_2048, sigma_LOS_2048), '--', 'Color', color_2048, 'LineWidth', 3, 'DisplayName', sprintf('Rician 2048 (K=%.1f dB)', K_dB_2048));
plot(r_rice, rician_pdf(r_rice, A_LOS_8192, sigma_LOS_8192), '--', 'Color', color_8192, 'LineWidth', 3, 'DisplayName', sprintf('Rician 8192 (K=%.1f dB)', K_dB_8192));

xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'LOS_Rician_Distribution.png'));

% Goodness-of-fit test: Kolmogorov-Smirnov test for Rician
fprintf('\n--- LOS Rician Distribution K-S Test ---\n');
[h_los_512, p_los_512, ks_los_512] = kstest(amp_LOS_512, 'CDF', makedist('Rician', 's', A_LOS_512, 'sigma', sigma_LOS_512));
fprintf('LOS 512:  K-S statistic = %.4f, p-value = %.4f %s\n', ks_los_512, p_los_512, fit_quality(ks_los_512));

[h_los_2048, p_los_2048, ks_los_2048] = kstest(amp_LOS_2048, 'CDF', makedist('Rician', 's', A_LOS_2048, 'sigma', sigma_LOS_2048));
fprintf('LOS 2048: K-S statistic = %.4f, p-value = %.4f %s\n', ks_los_2048, p_los_2048, fit_quality(ks_los_2048));

[h_los_8192, p_los_8192, ks_los_8192] = kstest(amp_LOS_8192, 'CDF', makedist('Rician', 's', A_LOS_8192, 'sigma', sigma_LOS_8192));
fprintf('LOS 8192: K-S statistic = %.4f, p-value = %.4f %s\n', ks_los_8192, p_los_8192, fit_quality(ks_los_8192));
fprintf('Note: K-S statistic interpretation: <0.05=excellent, <0.10=good, <0.20=reasonable\n');

%% -----------------------------------------------------------
% Individual Histogram Plots with PDF Overlay
% -----------------------------------------------------------

% NLOS 512 Histogram with Rayleigh PDF
figure(5);
clf;
hold on;
histogram(amp_NOLOS_512, 50, 'Normalization', 'pdf', 'FaceColor', color_512, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'NLOS 512');
r_ray = linspace(0, max(amp_NOLOS_512), 200);
plot(r_ray, rayleigh_pdf(r_ray, sigma_NOLOS_512), '-', 'Color', color_512, 'LineWidth', 3, 'DisplayName', 'Rayleigh PDF');
xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'NLOS_512_Histogram.png'));

% NLOS 2048 Histogram with Rayleigh PDF
figure(6);
clf;
hold on;
histogram(amp_NOLOS_2048, 50, 'Normalization', 'pdf', 'FaceColor', color_2048, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'NLOS 2048');
r_ray = linspace(0, max(amp_NOLOS_2048), 200);
plot(r_ray, rayleigh_pdf(r_ray, sigma_NOLOS_2048), '-', 'Color', color_2048, 'LineWidth', 3, 'DisplayName', 'Rayleigh PDF');
xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'NLOS_2048_Histogram.png'));

% NLOS 8192 Histogram with Rayleigh PDF
figure(7);
clf;
hold on;
histogram(amp_NOLOS_8192, 50, 'Normalization', 'pdf', 'FaceColor', color_8192, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'NLOS 8192');
r_ray = linspace(0, max(amp_NOLOS_8192), 200);
plot(r_ray, rayleigh_pdf(r_ray, sigma_NOLOS_8192), '-', 'Color', color_8192, 'LineWidth', 3, 'DisplayName', 'Rayleigh PDF');
xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'NLOS_8192_Histogram.png'));

% LOS 512 Histogram with Rician PDF
figure(8);
clf;
hold on;
histogram(amp_LOS_512, 50, 'Normalization', 'pdf', 'FaceColor', color_512, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'LOS 512');
r_rice = linspace(0, max(amp_LOS_512), 200);
plot(r_rice, rician_pdf(r_rice, A_LOS_512, sigma_LOS_512), '-', 'Color', color_512, 'LineWidth', 3, 'DisplayName', sprintf('Rician (K=%.1f dB)', K_dB_512));
xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'LOS_512_Histogram.png'));

% LOS 2048 Histogram with Rician PDF
figure(9);
clf;
hold on;
histogram(amp_LOS_2048, 50, 'Normalization', 'pdf', 'FaceColor', color_2048, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'LOS 2048');
r_rice = linspace(0, max(amp_LOS_2048), 200);
plot(r_rice, rician_pdf(r_rice, A_LOS_2048, sigma_LOS_2048), '-', 'Color', color_2048, 'LineWidth', 3, 'DisplayName', sprintf('Rician (K=%.1f dB)', K_dB_2048));
xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'LOS_2048_Histogram.png'));

% LOS 8192 Histogram with Rician PDF
figure(10);
clf;
hold on;
histogram(amp_LOS_8192, 50, 'Normalization', 'pdf', 'FaceColor', color_8192, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'LOS 8192');
r_rice = linspace(0, max(amp_LOS_8192), 200);
plot(r_rice, rician_pdf(r_rice, A_LOS_8192, sigma_LOS_8192), '-', 'Color', color_8192, 'LineWidth', 3, 'DisplayName', sprintf('Rician (K=%.1f dB)', K_dB_8192));
xlabel('Amplitude', 'FontSize', 32);
ylabel('PDF', 'FontSize', 32);
set(gca, 'FontSize', 32);
grid on;
legend('Location', 'best', 'FontSize', 32);
hold off;
saveas(gcf, fullfile(scriptDir, 'LOS_8192_Histogram.png'));

% Helper function to assess fit quality based on K-S statistic
function quality = fit_quality(ks_stat)
    if ks_stat < 0.05
        quality = '(excellent fit)';
    elseif ks_stat < 0.10
        quality = '(good fit)';
    elseif ks_stat < 0.20
        quality = '(reasonable fit)';
    else
        quality = '(poor fit)';
    end
end
