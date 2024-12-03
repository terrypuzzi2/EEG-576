clc;
clear all;

% Load the data
data1 = readtable('Green_1.csv'); % Replace with your first file name
data2 = readtable('Red_3.csv'); % Replace with your second file name

% Extract EEG signals
eegSignals1 = table2array(data1(4:end, 2:9)); % Columns 2-10: EEG signals from Green_1
eegSignals2 = table2array(data2(4:end, 2:9)); % Columns 2-10: EEG signals from Green_2
fs = 255; % Sampling frequency in Hz

% Define the time intervals
timeInterval1 = [15, 25]; % Time segment 1 in seconds
timeInterval2 = [35, 45]; % Time segment 2 in seconds

% Convert time intervals to sample indices
sampleInterval1 = round(timeInterval1 * fs);
sampleInterval2 = round(timeInterval2 * fs);

% Extract data for the specified time intervals
segment1_1 = eegSignals1(sampleInterval1(1):sampleInterval1(2), :);
segment1_2 = eegSignals2(sampleInterval1(1):sampleInterval1(2), :);

segment2_1 = eegSignals1(sampleInterval2(1):sampleInterval2(2), :);
segment2_2 = eegSignals2(sampleInterval2(1):sampleInterval2(2), :);

% Combine the two segments for both datasets
combinedData1 = [segment1_1; segment2_1];
combinedData2 = [segment1_2; segment2_2];

% Design a bandpass filter for alpha waves (8-13 Hz)
lowCutoff = 8; % Lower bound of alpha waves in Hz
highCutoff = 13; % Upper bound of alpha waves in Hz
[b, a] = butter(4, [lowCutoff, highCutoff] / (fs / 2), 'bandpass'); % 4th-order Butterworth filter

% Apply the bandpass filter to both datasets
filteredSignals1 = filtfilt(b, a, combinedData1);
filteredSignals2 = filtfilt(b, a, combinedData2);

% Calculate power spectral density (PSD) for each channel
[psd1, f] = pwelch(filteredSignals1, [], [], [], fs);
[psd2, ~] = pwelch(filteredSignals2, [], [], [], fs);

% Select the alpha wave range (8-13 Hz)
alphaRange = (f >= 8 & f <= 13);
psd1Alpha = psd1(alphaRange, :); % PSD values in alpha range for all channels
psd2Alpha = psd2(alphaRange, :); % PSD values in alpha range for all channels

% Initialize for correlation values
nChannels = size(filteredSignals1, 2);
correlationValues = zeros(1, nChannels); % To store correlation coefficients
maxLag = 10; % Define maximum lag in samples

% Loop through channels for cross-correlation and correlation coefficient
for channel = 1:nChannels
    % Compute cross-correlation for the alpha power of corresponding channels
    [xcorrResult, lags] = xcorr(psd1Alpha(:, channel), psd2Alpha(:, channel), maxLag, 'coeff');
    
    % Find the lag with the maximum correlation
    [maxCorr, idx] = max(xcorrResult);
    optimalLag = lags(idx) / fs; % Convert lag to seconds
    
    % % Display cross-correlation results
    % fprintf('Channel %d: Max Corr = %.2f at Lag = %.2f seconds\n', channel, maxCorr, optimalLag);

    % Compute general correlation coefficient
    R = corrcoef(psd1Alpha(:, channel), psd2Alpha(:, channel));
    correlationValues(channel) = R(1, 2); % Extract correlation coefficient
    
    % Display general correlation coefficient
    fprintf('Channel %d: General Correlation Coefficient = %.2f\n', channel, correlationValues(channel));
end
