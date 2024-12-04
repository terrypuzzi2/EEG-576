clc;
clear;
close all;

% Load the data
data = readtable('Green_1.csv'); % Replace with your file name

% Extract sample indices and EEG signals
sampleIndex = table2array(data(4:end, 1)); % Sample indices
eegSignals = table2array(data(4:end, 2:9)); % EEG signals, 8 channels
fs = 256; % Sampling frequency in Hz

% Define the time intervals
timeInterval1 = [15, 25]; % in seconds
timeInterval2 = [35, 45]; % in seconds

% Convert time intervals to sample indices
sampleInterval1 = round(timeInterval1 * fs);
sampleInterval2 = round(timeInterval2 * fs);

% Extract data for the specified time intervals
dataSegment1 = eegSignals(sampleInterval1(1):sampleInterval1(2), :);
dataSegment2 = eegSignals(sampleInterval2(1):sampleInterval2(2), :);

% Combine the two segments
combinedData = [dataSegment1; dataSegment2];

% Design a bandpass filter for alpha waves (8-13 Hz)
lowCutoff = 8; % Lower bound of alpha waves in Hz
highCutoff = 13; % Upper bound of alpha waves in Hz
[b, a] = butter(4, [lowCutoff, highCutoff] / (fs / 2), 'bandpass'); % 4th-order Butterworth filter

% Apply the filter to each channel
filteredSignals = filtfilt(b, a, combinedData);

% Calculate the Fourier Transform and Power Spectral Density (PSD)
n = size(filteredSignals, 1); % Number of samples in combined data
frequencies = (0:n-1) * (fs / n); % Frequency axis
alphaRange = (frequencies >= 8 & frequencies <= 13); % Alpha wave frequency range

% Initialize storage for PSD and frequency data
psdData = zeros(sum(alphaRange), size(filteredSignals, 2));

for channel = 1:size(filteredSignals, 2)
    % Perform FFT
    fftResult = fft(filteredSignals(:, channel));
    
    % Compute Power Spectrum
    powerSpectrum = abs(fftResult).^2 / n;
    
    % Extract the alpha range
    psdData(:, channel) = powerSpectrum(alphaRange);
end

% Plot the results
for channel = 1:size(filteredSignals, 2)
    figure;
    plot(frequencies(alphaRange), psdData(:, channel));
    title(['Channel ', num2str(channel), ' Alpha Waves (8-13 Hz)']);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    grid on;
end
