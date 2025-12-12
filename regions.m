% Region limits around the pyramidal horn antenna at 4.8 GHz
% Method per Chapter 6/IEEE: use D = largest physical dimension of aperture.

clear; clc;

% Parameters
f  = 4.8e9;      % Frequency in Hz
c  = 3e8;        % Speed of light in m/s
lambda = c / f;  % Wavelength in meters

% Aperture dimensions (meters)
a1 = 10.4e-2;    % Width  = 10.4 cm
b1 = 13.8e-2;    % Height = 13.8 cm

% Characteristic sizes
D = max(a1, b1);            % Largest physical dimension (used for region limits)
D_diag = sqrt(a1^2 + b1^2); % Diagonal (optional reference)

% Region boundaries (IEEE definitions)
% Reactive near-field upper limit: R >= 0.62 * sqrt(D^3 / lambda)
R_rnf = 0.62 * sqrt(D^3 / lambda); % R_fresnel

% Fresnel region: three criteria, take maximum
% 1. Wavelength-based: R > 20*λ/(2π)
% 2. Geometry-based: R > 20*ρ_max (ρ_max = distance from phase center to farthest aperture corner)
% 3. Rayleigh: R > 2D²/λ
rho_max = sqrt((a1/2)^2 + (b1/2)^2); % distance from center to farthest corner
R_wavelength = 20 * lambda / (2*pi);
R_geometry = 20 * rho_max;
R_fresnel = 0.62 * sqrt(D^3 / lambda);
R_fresnel_max = max([R_wavelength, R_geometry, R_fresnel]);
R_fraunhofer = 2 * D^2 / lambda;
R_fraunhofer_max = max([R_wavelength, R_geometry, R_fraunhofer]);
% Display
fprintf('Frequency               : %.2f GHz\n', f/1e9);
fprintf('Wavelength              : %.3f cm\n', lambda*100);
fprintf('Aperture (W x H)        : %.1f cm × %.1f cm\n', a1*100, b1*100);
fprintf('Largest side D          : %.1f cm\n', D*100);

fprintf('\n=== REGION BOUNDARIES ===\n');
fprintf('Reactive near-field     : R < %.1f cm\n', R_rnf*100);
fprintf('\nFresnel region:\n');
fprintf('  1. Wavelength (20λ/2π)  : R > %.1f cm\n', R_wavelength*100);
fprintf('  2. Geometry (20ρ_max)   : R > %.1f cm\n', R_geometry*100);
fprintf('  3. Fresnel (0.6 sqrt(D³/λ))     : R > %.1f cm\n', R_fresnel*100);
fprintf('\nFresnel boundary: R >= %.1f cm\n', R_fresnel_max*100);
fprintf('Far-field (Fraunhofer)  :\n');
fprintf('  1. Wavelength (20λ/2π)  : R > %.1f cm\n', R_wavelength*100);
fprintf('  2. Geometry (20ρ_max)   : R > %.1f cm\n', R_geometry*100);
fprintf('  3. Fraunhover (2D²/λ)     : R > %.1f cm\n', R_fraunhofer*100);
fprintf('\nFar-field boundary: R >= %.1f cm\n', R_fraunhofer_max*100);
% Field convergence check using R·|E| stabilization at multiple angles
% Set simulateFields = true to compute distance-normalized fields.
simulateFields = true;
if simulateFields
    freq = f;
    h = horn;
    h.FlareWidth = b1; h.FlareHeight = a1; h.FlareLength = 0.174;
    h.Width = 0.0475; h.Height = 0.022; h.Length = 0.05;
    
    
    % Distance sweep (meters)
    R_sweep = linspace(0.1, 3, 100); % 10 cm to 6 m
    
    % Angles from boresight (degrees) - in elevation plane
    angles_deg = [0, 10, 20, 30, 40];
    nAngles = numel(angles_deg);
    
    % Preallocate
    RE_normalized = zeros(nAngles, numel(R_sweep));
    
    fprintf('\n--- Simulated R·|E| Stabilization Analysis ---\n');
    
    % Compute R·|E| for each angle and distance
    % Horn boresight is along +x axis
    for iAng = 1:nAngles
        theta = angles_deg(iAng);
        for k = 1:numel(R_sweep)
            r = R_sweep(k);
            % Boresight +x, angle θ in xz-plane
            x = r * cosd(theta);  % along boresight
            z = r * sind(theta);  % off-axis
            p = [x; 0; z];
            [E, ~] = EHfields(h, freq, p);
            RE_normalized(iAng, k) = r * norm(E);
        end
    end
    
    % Plot R·|E| vs distance for each angle
    figure;
    hold on;
    colors = lines(nAngles);
    legendEntries = cell(nAngles, 1);
    
    for iAng = 1:nAngles
        plot(R_sweep*100, RE_normalized(iAng, :), '-o', 'Color', colors(iAng, :), 'LineWidth', 1.5);
        legendEntries{iAng} = sprintf('%d°', angles_deg(iAng));
    end
    
    % Add region boundaries
    yLimits = ylim;
    plot([R_fresnel*100 R_fresnel*100], yLimits, 'k-', 'LineWidth', 2);
    plot([R_fraunhofer*100 R_fraunhofer*100], yLimits, 'k-', 'LineWidth', 2);
    plot([R_fraunhofer_max*100 R_fraunhofer_max*100], yLimits, 'k-', 'LineWidth', 2);
    
    xlabel('Distance (cm)');
    ylabel('R·|E| (normalized field)');
    title('Far-field Convergence (MatLab Simulation): R·|E| vs Distance');
    legend(legendEntries, 'Location', 'best');
    grid on;
    hold off;
    
    % Find stabilization distance (within 1% of farthest value)
    tol = 0.01;
    fprintf('\nStabilization distances (R·|E| within 1%% of reference):\n');
    for iAng = 1:nAngles
        ref_val = RE_normalized(iAng, end);
        rel_err = abs(RE_normalized(iAng, :) - ref_val) / ref_val;
        idx = find(rel_err < tol, 1, 'first');
        if ~isempty(idx)
            fprintf('  %2d° off-axis: %.1f cm\n', angles_deg(iAng), R_sweep(idx)*100);
        else
            fprintf('  %2d° off-axis: not converged within sweep\n', angles_deg(iAng));
        end
    end
end
