%%  TUS-EEG processing. 

% SAFE-TUS: Signal Artifact Filtering for EEG with FUS


% DRCMR 2025 - XCT - BrainStim Methods group

% --> GENERAL INFO:
% The current example preprocessing pieplpile is intended to provide a simple and
% fast preprocessing pipeline focusing on removing DC-coupling like
% artifacts. Further artifacts (muscular, line noise, ocular...) may be
% considered.


% ----------------> STEPS:
% 0) Load data
% 1) Epoch
% 2) De-mean data
% 2.1) Remove bad trials and electrodes (optionally)
% 3) Remove TUS artifacts
% 4) Interpolate removed signal
% 5) Bandpass filtering 2-500Hz, Baseline correction and  Rereferencing to the average.
% 8) Save datasets

%% Addpaths

clear all;
close all;
clc

% Path to eeglab
addpath('.../eeglab2025.0.0/')

% Path to fieltrip
addpath('.../fieldtrip-20240110/')


%% 0. LOAD DATA

% Open EEGlab 
eeglab;

% Define file to load 
name_dataset = 'data.vhdr';
path_dataset = '.../eeg/';

% Load file
EEG = pop_loadbv(path_dataset, name_dataset);

eeglab redraw

%% 0.1 Preliminar processing Epoch the data  

% Remove non existing channels due to TUS transducer placement
EEG = pop_select( EEG, 'rmchannel',{'44','45','56','57','58'});
eeglab redraw

% Filter events for epochs
% 
 remove_idx = ismember({EEG.event.type}, {'S 15'});
 EEG.event = EEG.event(~remove_idx);
% 
 remove_idx = ismember({EEG.urevent.type}, {'S 15'});
 EEG.urevent = EEG.urevent(~remove_idx);
%
% Step 1: Filter events
keep_idx = ismember({EEG.event.type}, {'S 11'});
EEG.event = EEG.event(keep_idx);

% Step 2: Extract latencies
latencies = [EEG.event.latency];

% Step 3: Sort events by latency (in case they're not sorted)
[latencies, sort_idx] = sort(latencies);
EEG.event = EEG.event(sort_idx);

% Step 4: Compute latency differences to find stimulus trains
lat_diff = diff(latencies);
threshold = 1000; % <-- adjust this based on your actual ISI (inter-stimulus interval)
train_start_idx = [1, find(lat_diff > threshold) + 1];

% Step 5: Rename events
for i = 1:length(EEG.event)
    EEG.event(i).type = 'Stimuly'; % default to Stimuly
end

% Label first event of each train as 'S 11'
for idx = train_start_idx
    EEG.event(idx).type = 'S 11';
end

% Filter ur epochs
% Step 1: Filter urevents
keep_idx = ismember({EEG.urevent.type}, {'S 11'});
EEG.urevent = EEG.urevent(keep_idx);

% Step 2: Extract latencies
latencies = [EEG.urevent.latency];

% Step 3: Sort urevents by latency (in case they're not sorted)
[latencies, sort_idx] = sort(latencies);
EEG.urevent = EEG.urevent(sort_idx);

% Step 4: Compute latency differences to find stimulus trains
lat_diff = diff(latencies);
threshold = 1000; % <-- adjust this based on your actual ISI (inter-stimulus interval)
train_start_idx = [1, find(lat_diff > threshold) + 1];

% Step 5: Rename urevents
for i = 1:length(EEG.urevent)
    EEG.urevent(i).type = 'Stimuly'; % default to Stimuly
end

% Label first urevent of each train as 'S 11'
for idx = train_start_idx
    EEG.urevent(idx).type = 'S 11';
end

%eeglab redraw


%% 1. Epoch

EEG = pop_epoch( EEG, {  'S 11'  }, [-2        2]); % Pulses are stamped on the EEG as 'S 15' markers. Modify your marker if necessary for epoching.
%eeglab redraw

%Figure
figure;
plot(EEG.times, EEG.data(:,:,4), 'b'); 
%xlim([-100 2000]); 
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% 2. Demean 

EEG = pop_rmbase( EEG, [-2000 1999]); %selet window according to Epoch length


% Figure
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-100 2000]); 
%ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


% 2.1 RENAME TRIGGERS TO REMOVE PULSES
%  Rename events
for i = 1:length(EEG.event)
    EEG.event(i).type = 'S 15'; % default to Stimuly
end

% Rename urevents
for i = 1:length(EEG.urevent)
    EEG.urevent(i).type = 'S 15'; % default to Stimuly
end

%% Remove bad channels manually if necessary 

% EEG = pop_select( EEG, 'rmchannel',{'46'});



%% Remove bad trials
%% Trials rejection -(identify bad trials) --> Select the bad trials

TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'
 
%% Trials rejection - (remove bad trials)

if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

eeglab redraw


%% 3. Interpolate sharp TUS pulses given  hardware (ACTICHAMP) filtering of DC artifacts. Withoud DC filering, artifacts may look like DC sharp step like transitions.

% Code assuming TUS pulse 10ms
EEG = pop_tesa_removedata( EEG, [-1.5,1.5], [], {'S 15'} ); %remove from -1.5 to 1.5 ms to remove the up-phase of the DC
EEG = pop_tesa_removedata( EEG, [3.5, 6.5], [], {'S 15'} ); %remove from -ms to remove the up-phase of the down. Adjust according to TUS PRF

% eeglab redraw

% Plot data
figure; pop_timtopo(EEG, [-20 200], [-20 -5 3 9 99 103 106 110 ], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-7 12]);  ylim([-13 13])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 4. INTERPOLATE DATA

EEG = pop_tesa_interpdata( EEG, 'cubic', [50,50] );

%eeglab redraw

% Plot data
figure; pop_timtopo(EEG, [-2000 1999], [-250 -500 500 150 250], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-2000 2000]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


% Filter epochs again
% Step 1: Filter events
keep_idx = ismember({EEG.event.type}, {'S 11', 'S 15'});
EEG.event = EEG.event(keep_idx);

% Step 2: Extract latencies
latencies = [EEG.event.latency];

% Step 3: Sort events by latency (in case they're not sorted)
[latencies, sort_idx] = sort(latencies);
EEG.event = EEG.event(sort_idx);

% Step 4: Compute latency differences to find stimulus trains
lat_diff = diff(latencies);
threshold = 1000; % <-- adjust this based on your actual ISI (inter-stimulus interval)
train_start_idx = [1, find(lat_diff > threshold) + 1];

% Step 5: Rename events
for i = 1:length(EEG.event)
    EEG.event(i).type = 'Stimuly'; % default to Stimuly
end

% Label first event of each train as 'S 11'
for idx = train_start_idx
    EEG.event(idx).type = 'S 11';
end



%% 5.  Filtering, baseline correction and re-referencing to the average

%fitler
 EEG = pop_tesa_filtbutter(EEG, 2, [], 4, 'highpass');
 EEG = tesa_filtbutter(EEG, [], 500, 4, 'lowpass');

%baseline correction
EEG = pop_rmbase(EEG, [-210 -10]); %ms

%re-reference
EEG = pop_reref(EEG, []);

eeglab redraw



% Plot data
figure; pop_timtopo(EEG, [-2000 1999], [-250 -500 500 150 250], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-2000 2000]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');



%% 6. SAVE DATASET

EEG = pop_saveset( EEG, 'filename',[name_dataset,'_clean.set'],'filepath',[path_dataset]);

%% END

