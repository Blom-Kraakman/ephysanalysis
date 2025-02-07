%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   align spikes to event timing

clearvars

% TODO: parse out mouse number and recording site to avoid mistakes when
% setting directories

% set directories
recordingFolder = 'D:\DATA\EphysRecordings\M11\M11_2024-05-22_12-34-22\';
recPath = [recordingFolder 'Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\'];
TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
KSPath = 'D:\DATA\EphysRecordingsSorted\M11\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M11\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M11'; % output directory

rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording

Fs = 30000; % sampling freq

relevant_sessions = [1 10];
skip_sessions = [];

%relevant_sessions = [7 8]; %M12 ICX 1:4, 5:9; ICC 10:13
%relevant_sessions = [10 13]; %M13 ICX 1:6; ICC 10:13
%skip_sessions = [1 2 3 4 5 6 9]; %M13.1 7:9
%relevant_sessions = [1 11]; % M8
%skip_sessions = 10; % M8
%relevant_sessions = [1 10]; %M11 and M10
%skip_sessions = 2;
%relevant_sessions = [1 9]; % M6 behaviour files (if only 1 behavior file in rec: [1 1])
%relevant_sessions = [4 11]; % M9

% save session TTLs if needed, else load correct file
% get trials onset TTLs of all sessions in recording
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']);

if isempty(TTLs_file)
    % get trials onset TTLs of all sessions in recording
    [sessions_TTLs, sessions_TTLs_variables] = getSessionTTLs(messagesPath, rec_samples, Fs, skip_sessions);
    
    % save
    filename = sprintf('M%s_S%02d-%02d_OE_TTLs', BehaviorPath(29:30), relevant_sessions(1), relevant_sessions(2));
    save(fullfile(OutPath, filename), "sessions_TTLs", "sessions_TTLs_variables")
    disp('sessions TTLs saved')
else % load from directory
    TTLs = load([OutPath '\' TTLs_file.name]); 
    disp('sessions TTLs loaded from saved file')
    sessions_TTLs = TTLs.sessions_TTLs;
end

%% Kilosort: post-curation unit extraction
% IronClust: post-curation unit extraction [spiketimes, cids,cpos] = ircGoodClusters(spiketimecsv,clusterqualitycsv);
% spike extraction from curated units

% to add: load if previously saved
[spiketimes, cids] = extractspikes(BehaviorPath, KSPath, TTLPath, relevant_sessions, rec_samples, Fs, OutPath);

% align spikes
% saves aligned spikes with corresponding stimulus onset and offset timings
% spike times in sec, Srise & Sfall in samples

alignspikes(BehaviorPath, TTLPath, OutPath, spiketimes, relevant_sessions, skip_sessions, cids, sessions_TTLs, Fs);

% ----------------------- FRA analysis & Plotting ----------------------- %%
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster

close all

% load FRA session(s) aligned spikes
aligned_spikes_files = dir(fullfile(OutPath, '*FRA_AlignedSpikes.mat'));

for file = 1:size(aligned_spikes_files, 1)

    % load each FRA aligned spikes file
    aligned_spikes = load([aligned_spikes_files(file).folder '\' aligned_spikes_files(file).name]);

    % load corresponding stimuli file
    session = aligned_spikes_files(file).name(6:7);
    stim_files = dir(fullfile(BehaviorPath, ['*_S' session '_*.mat']));

    if size(stim_files, 1) ~= 1 % check input
        error('Mismatch between stimulus and aligned spikes file')
    end

    stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

    % FRA analysis saves heatmap figures
    FSL = 0;
    FRAanalysis(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, FSL);

end

%% optional: match units between session
% 
% OutPath = 'D:\DATA\Processed\M16'; % output directory
% 
% rec1 = [273 269 287 279 201 298 212 318 321 238 256 258];
% rec2 = [151 153 45 53 69 78 237 163 184 196 175 201];
% matchedUnits = [rec1', rec2'];
% 
% % save
% filename = sprintf('M16_ICX_MatchedUnits');
% save(fullfile(OutPath, filename), "matchedUnits")
