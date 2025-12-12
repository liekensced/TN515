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
set(groot, 'DefaultAxesFontSize', 18); % bump plot font sizes
set(groot, 'DefaultTextFontSize', 20);


% Set the geometric properties
h = horn;
h.FlareWidth = 0.138; h.FlareHeight = 0.104; h.FlareLength = 0.174;
h.Width = 0.0475; h.Height = 0.022; h.Length = 0.05;
freq = 4.8e9;  % 4.8 GHz


% Display antenna info
figure;
show(h);
title('Pyramidal Horn Antenna Geometry');


% % Plot 3D radiation pattern
% figure;
% pattern(h, freq);
% [pat, az, el] = pattern(h, freq);
% title('3D Radiation Pattern at 4.8 GHz');
% Gpeak = max(pat(:));

% % Plot azimuth (H-plane) pattern
% figure;
% patternAzimuth(h, freq, 0);  % Phi=0 degrees
% title('Azimuth (H-plane) Pattern at 4.8 GHz');

% % Plot elevation (E-plane) pattern
% figure;
% patternElevation(h, freq, 90);  % Phi=90 degrees
% title('Elevation (E-plane) Pattern at 4.8 GHz');


% % Far-field approximation distance
% D = max(h.FlareWidth, h.FlareHeight);  % Max dimension
% lambda = physconst('LightSpeed') / freq;
% R_far = 2 * D^2 / lambda;
% fprintf('Recommended far-field distance: %.2f m\n', R_far);