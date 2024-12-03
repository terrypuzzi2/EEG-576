clc
clear
close all

%% Root directory and file names
% root = 'your_directory_path'; % Replace with your directory path
files = {
    'BrainFlow-RAW_2024-11-14_13-17-10_4.csv',
    'BrainFlow-RAW_2024-11-14_13-17-10_5.csv',
    'BrainFlow-RAW_2024-11-14_13-17-10_6.csv'
};

% Sampling frequency
freq = 256; % Sampling frequency

% Define EEG bands and names
bands = {
    'Delta', [0.5, 4];   % Delta band
    'Theta', [4, 7];     % Theta band
    'Alpha', [8, 12];    % Alpha band
    'Beta', [13, 30];    % Beta band
    'Gamma', [30, 80]    % Gamma band
};

% Loop through each dataset
for i = 1:length(files)
    % Load the dataset
    data = load(fullfile(root, files{i}));
    
    % Extract sample numbers and EEG data
    sample_numbers = data(:, 1); % Time data
    eeg_data = data(:, 2:9);     % EEG data (columns 2 to 9)
    
    %% Remove faulty data from the start
    start_index = find(sample_numbers == 0, 1, 'first');
    sample_numbers = sample_numbers(start_index:end);
    eeg_data = eeg_data(start_index:end, :);

    %% Calculate duration of sample
    total_samples = length(sample_numbers);
    total_seconds = total_samples / freq; % Total duration in seconds

    % Generate continuous time vector
    time = (sample_numbers + floor((0:(total_samples - 1)) / freq) * freq) / freq;

    %% Limit the EEG data to a maximum of 55 seconds
    max_seconds = 55;
    max_samples = min(55 * freq, total_samples);
    time_limited = time(1:max_samples);
    eeg_data_limited = eeg_data(1:max_samples, :);

    %% Preprocessing for each band
    for b = 1:size(bands, 1)
        band_name = bands{b, 1};
        band_range = bands{b, 2};

        % Bandpass filter for the current band
        [b_coeff, a_coeff] = butter(2, band_range / (freq / 2), 'bandpass');
        filtered_EEG = filtfilt(b_coeff, a_coeff, eeg_data_limited);

        %% Save EEG data for EEGLAB
        % Transpose filtered_EEG to have channels as rows
        filtered_EEG = filtered_EEG';

        % Initialize EEGLAB dataset
        EEG = pop_importdata('dataformat', 'array', 'data', filtered_EEG, 'srate', freq, 'xmin', 0);

        % Set the number of channels and their labels
        EEG.nbchan = size(filtered_EEG, 1); % Number of rows in filtered_EEG
        EEG.chanlocs = struct('labels', {'Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5', 'Channel 6', 'Channel 7', 'Channel 8'});

        % Set output directory and dataset name
        EEG.filepath = pwd; % Current working directory
        EEG.setname = sprintf('%s_eeg_data_%d', band_name, i); % Unique name for each dataset

        % Save the dataset
        EEG = pop_saveset(EEG, 'filename', sprintf('%s_eeg_data_%d.set', band_name, i));

        %% Perform ICA decomposition
        EEG = pop_runica(EEG, 'extended', 1); % Use extended ICA
        EEG = pop_saveset(EEG, 'filename', sprintf('%s_eeg_ica_%d.set', band_name, i));

        %% Automatic artifact rejection
        EEG = pop_iclabel(EEG, 'default');
        EEG = pop_icflag(EEG, [NaN, NaN, 0.8, NaN, NaN, NaN, NaN]); % Retain non-artifact ICs

        % Save cleaned data
        pop_saveset(EEG, 'filename', sprintf('%s_eeg_cleaned_%d.set', band_name, i));
    end
end
