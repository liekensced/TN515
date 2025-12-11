% Quick plot of distance vs received power (auto-detects column names)
opts = detectImportOptions('measurements.csv');
% Use first two columns as distance and rx; detectImportOptions sanitizes headers
distVar = opts.VariableNames{1};
rxVar   = opts.VariableNames{2};
opts = setvartype(opts, {distVar, rxVar}, 'double');
T = readtable('measurements.csv', opts);

% Remove rows with missing distance or rx
valid = ~isnan(T.(distVar)) & ~isnan(T.(rxVar));
T = T(valid, :);

figure;
plot(T.(distVar), T.(rxVar), '-o', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Rx (dBm)');
grid on;
title('Free-space Rx vs Distance');
set(gca, 'XDir', 'reverse'); % optional: distances recorded from far to near
