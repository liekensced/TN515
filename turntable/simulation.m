%% Horn Antenna Simulation at 4.8 GHz
% This script simulates a pyramidal horn antenna using MATLAB Antenna Toolbox
% Dimensions (converted to meters):
% - Aperture width (FlareWidth): 10.4 cm = 0.104 m
% - Aperture height (FlareHeight): 13.8 cm = 0.138 m
% - Inner rib: 18.5 cm = 0.185 m -> Flare length = 0.174 m
% - Throat/waveguide width (Width): 2.2 cm = 0.022 m
% - Throat/waveguide height (Height): 4.75 cm = 0.0475 m
% Frequency: 4.8 GHz

clear; clc; close all;

% Create default horn antenna
h = horn;

% Set the geometric properties
h.FlareWidth = 0.104;   % Aperture width (m)
h.FlareHeight = 0.138;  % Aperture height (m)
h.FlareLength = 0.174;  % Flare length (m)
h.Width = 0.022;        % Waveguide/throat width (m)
h.Height = 0.0475;      % Waveguide/throat height (m)

% Other defaults (can adjust if needed)
% h.Length = 0.09;      % Waveguide length (default or set if known)
% h.FeedOffset = [0 0]; % Feed offset (default)

% Design the horn for 4.8 GHz operation
freq = 4.8e9;  % 4.8 GHz
h = design(h, freq);

% Display antenna info
%figure;
%show(h);
%title('Pyramidal Horn Antenna Geometry');


% Plot 3D radiation pattern
figure;
pattern(h, freq);
[pat, az, el] = pattern(h, freq);
title('3D Radiation Pattern at 4.8 GHz');
Gpeak = max(pat(:));

% Plot azimuth (H-plane) pattern
figure;
patternAzimuth(h, freq, 0);  % Phi=0 degrees
title('Azimuth (H-plane) Pattern at 4.8 GHz');

% Plot elevation (E-plane) pattern
figure;
patternElevation(h, freq, 90);  % Phi=90 degrees
title('Elevation (E-plane) Pattern at 4.8 GHz');


% Far-field approximation distance
D = max(h.FlareWidth, h.FlareHeight);  % Max dimension
lambda = physconst('LightSpeed') / freq;
R_far = 2 * D^2 / lambda;
fprintf('Recommended far-field distance: %.2f m\n', R_far);