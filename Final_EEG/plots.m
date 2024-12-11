clc;
clear;
close all;

newDataFile = 'Green_3_4.csv'; 
channelsToCorrelate = [2, 4];
newtimeInterval = [15, 25];

%% Set parameters
fs = 256; % Sampling frequency in Hz
lowCutoff = 32; % Lower bound in Hz
highCutoff = 100; % Upper bound in Hz

%% File names for the existing datasets
fileNames = { 'Green_2_1.csv', 'Green_2_2.csv','Green_2_3.csv', ...
              'Red_2_1.csv', 'Red_2_2.csv', 'Red_2_3.csv'};

%% Time intervals
timeIntervals = {[15, 25], [35, 45]}; % Two time intervals in seconds

%% Colors for the plots
colorMap = struct('Green', 'g', 'Red', 'r');

%% Initialize variables for storing PSD data for ANOVA
psdValuesChannels = cell(1, 6); % Cell array to store PSD values for each channel
categoriesChannels = cell(1, 6); % Cell array to store categories for each channel
processedData = struct();

%% ICA File Check for Existing Data
icaFileName = 'icaData.mat'; % File to store ICA results

% Check if ICA data already exists
if exist(icaFileName, 'file')
    disp('Loading ICA data from file...');
    load(icaFileName, 'icaData'); % Load precomputed ICA data
else
    disp('Performing ICA and saving results...');
    icaData = struct(); % Initialize a structure to hold ICA data
    % Loop through each channel (1 to 6)
    for channel = 1:6
        for i = 1:length(fileNames)
            % Load the data
            data = readtable(fileNames{i});
            eegSignalsPreICA = table2array(data(4:end, [2:3, 6:9])); % EEG signals, 6 channels (skip 4 and 5)
            

            %% Perform ICA
            [weights, sphere, a, b, c, d] = runica(eegSignalsPreICA', 'pca', 6); % Assuming 6 components
            icaSignals = weights * eegSignalsPreICA'; % ICA decomposition, transposed for proper matrix multiplication
            
            % Store ICA data with more descriptive labels
            icaData(i).(sprintf('Channel_%d', channel)).signals = icaSignals'; % Store ICA signals for each file/channel
            icaData(i).(sprintf('Channel_%d', channel)).fileName = fileNames{i}; % Store corresponding file name
            icaData(i).(sprintf('Channel_%d', channel)).channel = channel; % Store channel number
        end
    end
    % Save the ICA data to file for future use
    save(icaFileName, 'icaData');
end

%% Process the data and compute PSD for each file using ICA results
for channel = 1:6
    for i = 1:length(fileNames)
        % Retrieve the ICA data for the current file and channel
        icaSignals = icaData(i).(sprintf('Channel_%d', channel)).signals;

        %% Determine the category and color based on the filename
        if contains(fileNames{i}, 'Green')
            category = 'Green';
            plotColor = colorMap.Green;
        elseif contains(fileNames{i}, 'Red')
            category = 'Red';
            plotColor = colorMap.Red;
        end
        
        %% Process each time interval
        for j = 1:length(timeIntervals)
            %% Convert time interval to sample indices
            timeInterval = timeIntervals{j};
            sampleInterval = round(timeInterval * fs);

            %% Extract data for the specified time interval
            dataSegment = icaSignals(sampleInterval(1):sampleInterval(2), :);

            %% Design a bandpass filter
            [b, a] = butter(4, [lowCutoff, highCutoff] / (fs / 2), 'bandpass');

            %% Apply the filter to the current channel
            filteredSignal = filtfilt(b, a, dataSegment(:, channel));

            %% Fourier Transform and Power Spectral Density (PSD)
            n = size(filteredSignal, 1); % Number of samples in the segment
            frequencies = (0:n-1) * (fs / n); % Frequency axis
            alphaRange = (frequencies >= lowCutoff & frequencies <= highCutoff); % Alpha range

            %% Calculate Power Spectrum for the channel
            fftResult = fft(filteredSignal);
            powerSpectrum = abs(fftResult).^2 / n; % Power spectrum

            %% Extract the alpha range and store PSD
            alphaPSD = powerSpectrum(alphaRange);
            psdValuesChannels{channel} = [psdValuesChannels{channel}; mean(alphaPSD)];
            categoriesChannels{channel} = [categoriesChannels{channel}; {category}];

            %% Replace invalid characters in the key
            key = sprintf('%s_Interval%d_Channel%d', strrep(fileNames{i}, '.', '_'), j, channel);
            processedData.(key).PSD = alphaPSD;
            processedData.(key).Category = category;
        end
    end
end

%% Perform ANOVA for each channel
for channel = 1:6
    [p, tbl, stats] = anova1(psdValuesChannels{channel}, categoriesChannels{channel}, 'off');
    
    if p < 0.05
        fprintf('\nANOVA Results for Channel %d:\n', channel);
        disp(['p-value: ', num2str(p)]);
        disp('Significant difference between Red and Green categories.');
    end
end

%% --- New Dataset for Prediction ---
newData = readtable(newDataFile);
newEEGSignalsPreICA = table2array(newData(4:end, [2:3, 6:9])); % EEG signals, 6 channels

% Perform ICA decomposition on the new data (testing dataset)
[weights, sphere, ~, ~, ~, ~] = runica(newEEGSignalsPreICA', 'pca', 6); % Assuming 6 components
icaSignalsNew = weights * newEEGSignalsPreICA'; % ICA decomposition

% Convert ICA signals back to EEG data (for visualization or further processing)
newEEGSignals = icaSignalsNew';

%% Process channels for the 15-17 second interval
results = struct();

for channelIdx = 1:length(channelsToCorrelate)
    channel = channelsToCorrelate(channelIdx);
    sampleInterval = round(newtimeInterval * fs);

    %% Extract data for the interval
    newDataSegment = newEEGSignals(sampleInterval(1):sampleInterval(2), channel);

    %% Apply the same bandpass filter
    [b, a] = butter(4, [lowCutoff, highCutoff] / (fs / 2), 'bandpass');
    filteredNewData = filtfilt(b, a, newDataSegment);

    %% Fourier Transform and Power Spectral Density (PSD)
    n = size(filteredNewData, 1);
    frequencies = (0:n-1) * (fs / n);
    alphaRange = (frequencies >= lowCutoff & frequencies <= highCutoff);

    fftResult = fft(filteredNewData);
    powerSpectrum = abs(fftResult).^2 / n;
    newAlphaPSD = powerSpectrum(alphaRange);

    %% Correlate the selected channel's PSD to previous datasets
    correlationResults = struct();
    for key = fieldnames(processedData)'
        key = key{1};
        if contains(key, sprintf('Channel%d', channel))
            existingPSD = processedData.(key).PSD;

            % Directly calculate correlation (no need to align lengths)
            correlation = corr(newAlphaPSD, existingPSD, 'type', 'Kendall');
            correlationResults.(key) = correlation;
        end
    end

    %% Categorize correlations and calculate averages
    greenCorrelationSum = 0;
    greenCount = 0;
    redCorrelationSum = 0;
    redCount = 0;

    %% Get the keys from correlationResults
    correlationKeys = fieldnames(correlationResults);

    %% Ensure correlationResults is not empty
    if isempty(correlationKeys)
        error('No correlations were found. Check the processing of the PSD data and ensure keys exist.');
    end

    %% Loop through the keys
    for k = 1:length(correlationKeys)
        key = correlationKeys{k}; % Access key
        correlationValue = correlationResults.(key); % Access correlation value

        %% Categorize correlations based on key
        if contains(key, 'Green')
            greenCorrelationSum = greenCorrelationSum + correlationValue;
            greenCount = greenCount + 1;
        elseif contains(key, 'Red')
            redCorrelationSum = redCorrelationSum + correlationValue;
            redCount = redCount + 1;
        end
    end

    %% Calculate average correlations
    avgCorrelationGreen = greenCorrelationSum / max(greenCount, 1);
    avgCorrelationRed = redCorrelationSum / max(redCount, 1);

    %% Make prediction
    if avgCorrelationGreen > avgCorrelationRed
        predictedCategory = 'Green';
    else
        predictedCategory = 'Red';
    end

    %% Store results for this channel
    results.(sprintf('Channel%d', channel)) = predictedCategory;
end

%% Display prediction results
disp('Prediction results for the new dataset:');
disp(results);
