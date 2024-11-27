clc
clear
close all

%% Get data into matlab
directory=pwd;
root=directory;


data = load ('BrainFlow-RAW_2024-11-14_13-17-10_6.csv');
sample_numbers = data(:,1); %time data
eeg_data = data(:,2:9); %eeg data

%% Remove the faulty data from the start
start_index = find(sample_numbers == 0,1, 'first');
fsample_numbers = sample_numbers(start_index:end);
eeg_data = eeg_data(start_index:end, :);

%% Calculate duration of sample

freq = 256; %Sampling frequency
total_samples = length(sample_numbers);
total_seconds = total_samples/freq; %Total duration of the file in seconds

%Generate continuous time vector
time = (sample_numbers+floor((0:(total_samples - 1))/freq)*freq)/freq;




%% Minimizing the EEG data to max of 55 seconds

%Limit data for plotting
max_seconds = 55;
max_samples = min(55*freq,total_samples);
time_limited = time(1:max_samples);
eeg_data_limited = eeg_data(1:max_samples, :);


%% Preprocessing

%filtering for specific wave type (theta)
low_cutoff = 4 %low cutoff (Hz)
high_cutoff = 8 %high cutoff (Hz)
[b,a] = butter(2, [low_cutoff, high_cutoff]/(freq/2), 'bandpass');
filtered_EEG = filtfilt (b, a, eeg_data_limited);

%% Save EEG for EEGLAB

% Transpose filtered_EEG so that channels are rows
filtered_EEG = filtered_EEG';

% Correctly initialize the EEGLAB dataset
EEG = pop_importdata('dataformat', 'array', 'data', filtered_EEG, 'srate', freq, 'xmin', 0);

% Set the number of channels and channel labels
EEG.nbchan = size(filtered_EEG, 1); % Number of rows in filtered_EEG
EEG.chanlocs = struct('labels', {'Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5', 'Channel 6', 'Channel 7', 'Channel 8'});

% Specify a valid, writable directory
EEG.filepath = pwd; % Current working directory, or replace with a specific directory
EEG.setname = 'eeg_data'; % Name of the dataset

% Save the dataset
EEG = pop_saveset(EEG, 'filename', 'eeg_data.set');

%% EEGlab

%eeglab;

%ICA decomposition (filter out artifacts)

EEG= pop_runica(EEG, 'extended', 1); %use extended ICA
EEG= pop_saveset(EEG,'filename', 'eeg_ica.set');

EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG, [NaN, NaN, 0.8, NaN, NaN, NaN, NaN]); % Retain only non-artifact ICs

%%How to use EEGLAB GUI for artifact rejection
% 1. Go to Tools > Reject components by map.
% 2. Inspect each independent component's scalp map, power spectrum, and time series.
% 3. Identify and reject artifact-related components (e.g., eye blinks, muscle noise).
% 4. Save the cleaned dataset.

%Save data
pop_saveset(EEG, 'filename','eeg_cleaned.set');

%% Plotting the graph

%figure();
%hold on
%plot(time_limited, filtered_EEG);
%xlabel('Time (s)');
%xlim([0,max(time_limited)]);
%ylabel('EEG Amplitude (\muV)');
%title('EEG Data Over Time');
%legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5', 'Channel 6', 'Channel 7', 'Channel 8');
