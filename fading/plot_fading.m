% Plot fading dBm samples produced by GNU Radio
% Assumes 32-bit float little-endian samples (dBm) written raw to file.

fileName = 'LOS_2048_Sam.txt'; % change to any of the *_Sam files here
sampleRateHz = 1.28e6;           % capture sample speed (Hz)
fftSize = 8192;                  % FFT size used in GNU Radio

bytesPerSample = 4;              % float32
info = dir(fileName);
if isempty(info)
    error('File not found: %s', fileName);
end
if mod(info.bytes, bytesPerSample) ~= 0
    error('File size (%d bytes) is not divisible by %d; adjust format.', info.bytes, bytesPerSample);
end

fid = fopen(fileName, 'rb');
if fid < 0
    error('Could not open file: %s', fileName);
end
% Read float32 little-endian; change ''l'' to ''b'' if data are big-endian
samples = fread(fid, inf, 'float32', 0, 'l');
fclose(fid);

frameRate = sampleRateHz / fftSize;  % frames per second (power updates)
numSamples = numel(samples);
timeVec = (0:numSamples-1).' / frameRate;

titleStr = sprintf('%s (%d samples, %.2f s)', fileName, numSamples, timeVec(end));

figure;
plot(timeVec, samples, '-');
xlabel('Time (s)');
ylabel('Rx (dBm)');
title(titleStr);
grid on;

% Optional: quick smoothing for visualization (uncomment to use)
% win = 5; % samples
% hold on;
% plot(timeVec, movmean(samples, win), 'r', 'LineWidth', 1.2);
% legend('Raw', sprintf('Moving mean (%d)', win), 'Location', 'best');
