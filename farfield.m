% Far-field distance calculation for a pyramidal horn antenna
% Frequency: 4.8 GHz
% Aperture: 10.4 cm × 13.8 cm

clear; clc;

% Parameters
f  = 4.8e9;              % Frequency in Hz
c  = 3e8;                % Speed of light in m/s
lambda = c / f;          % Wavelength in meters

% Aperture dimensions in meters
a1 = 10.4e-2;            % Width  = 10.4 cm
b1 = 13.8e-2;            % Height = 13.8 cm  ← fixed line

% Largest dimension D (most accurate = diagonal of the aperture)
D = sqrt(a1^2 + b1^2);

% Far-field distance (standard formula)
R_far = 2 * D^2 / lambda;

% Display results
fprintf('Frequency         : %.2f GHz\n', f/1e9);
fprintf('Wavelength        : %.3f cm\n', lambda*100);
fprintf('Aperture          : %.1f cm × %.1f cm\n', a1*100, b1*100);
fprintf('Aperture diagonal : %.2f cm\n', D*100);
fprintf('\nFar-field distance (2D²/λ):\n');
fprintf('   %.1f cm  (%.2f m)\n', R_far*100, R_far);

% Optional: also show the value if only the longer side is used
R_far_longer_side = 2 * b1^2 / lambda;
fprintf('\nIf only the longer side (13.8 cm) is used: %.1f cm\n', R_far_longer_side*100);