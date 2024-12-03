%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   align spikes to event timing

clearvars

% set directories
recordingFolder = 'D:\DATA\EphysRecordings\M16\M16_2024-09-17_13-20-35\';
recPath = [recordingFolder 'Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\'];
TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
KSPath = 'D:\DATA\EphysRecordingsSorted\M16\ICX\rec2\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M16\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M16\ICX\rec2'; % output directory

rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording

Fs = 30000; % sampling freq

relevant_sessions = [6 12];
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

%% sessions TTLs
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

%% align spikes
% saves aligned spikes with corresponding stimulus onset and offset timings
% spike times in sec, Srise & Sfall in samples

alignspikes(BehaviorPath, TTLPath, OutPath, spiketimes, relevant_sessions, skip_sessions, cids, sessions_TTLs, Fs);

%% optional: match units between session

OutPath = 'D:\DATA\Processed\M16'; % output directory

rec1 = [273 269 287 279 201 298 212 318 321 238 256 258];
rec2 = [151 153 45 53 69 78 237 163 184 196 175 201];
matchedUnits = [rec1', rec2'];

% save
filename = sprintf('M16_ICX_MatchedUnits');
save(fullfile(OutPath, filename), "matchedUnits")

%% ----------------------- FRA analysis & Plotting ----------------------- %%

%% FRA analysis
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster
close all
OutPath = 'D:\DATA\Processed\M15\ICC';
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M15\';

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

% select correct input files
aligned_spikes = load([OutPath, '\M15_S16_FRA_AlignedSpikes']);
stimuli_parameters = load([BehaviorPath 'M15_S16_FRA.mat']);

% function saves figures
FSL = 0;
FRAanalysis(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, FSL);
