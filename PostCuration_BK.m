%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   align spikes to event timing

clearvars

% TODO: parse out mouse number and recording site to avoid mistakes when
% setting directories

% % set directories
% recordingFolder = 'D:\DATA\EphysRecordings\M27\M27_2025-05-27_12-26-52\';
% recPath = [recordingFolder 'Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\'];
% TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
% messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
% KSPath = 'D:\DATA\EphysRecordingsSorted\M27\'; % kilosort ephys data
% BehaviorPath = 'D:\DATA\Behavioral Stimuli\M27\'; % stimuli parameters
% OutPath = 'D:\DATA\Processed\M27'; % output directory

% set directories
recordingFolder = '\\store\department\neuw\shared\Aaron Wong\Data\InVivoEphys\Blom\EphysRecordings\M19\M19_2024-12-17_11-39-47_ICX\';
recPath = [recordingFolder 'Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\'];
TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
KSPath = '\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Blom\EphysRecordingsSorted\M19\ICX\'; % kilosort ephys data
BehaviorPath = '\\store\department\neuw\shared\Aaron Wong\Data\InVivoEphys\Blom\BehavioralStimuli\M19\'; % stimuli parameters
OutPath = '\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Blom\Processed\M19\ICX_test'; % output directory

Fs = 30000; % sampling freq

relevant_sessions = [1 17];
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
    rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording. supressed if alignedspikes done to save time
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

%  --- NEW ---
% loads data if previously saved
% saves cids and spike times in seperate files
spiketimes_file = dir([OutPath '\*_SpikeTimes.mat']);
if isempty(spiketimes_file)
    rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording. supressed if alignedspikes done to save time
    [spiketimes, ~] = extractspikes(BehaviorPath, KSPath, TTLPath, relevant_sessions, rec_samples, Fs, OutPath);
    disp('spiketimes and single cluster details saved')
else
    spiketimesdata = load([OutPath '\' spiketimes_file.name]);
    spiketimes = spiketimesdata.spiketimes;
    disp('spiketimes.mat loaded from saved files')
end

%% align spikes
% saves aligned spikes with corresponding stimulus onset and offset timings
% spike times in sec, Srise & Sfall in samples
% NEW: seperated function to call when aligning data

for session = relevant_sessions(1):relevant_sessions(2)

    % load data
    [~, stimuli_parameters, aligned_spikes, ~, ~, sessions_TTLs, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

    % skip earlier stim files not in rec session
    if isempty(dir([OutPath '\*_' sprintf('S%.2d_' , 4) '*_AlignedSpikes.mat' ]))
        disp(['Aligning session: ' num2str(session)])
    elseif ismember(str2double(stimuli_parameters.Par.Set), skip_sessions) || ~ismember(str2double(stimuli_parameters.Par.Set), relevant_sessions(1):relevant_sessions(2))
        continue
    else
        continue
    end

    % NEW: seperated generation of Srise from spike alignment to control moment to align to
    tsessions_TTLs = sessions_TTLs(sessions_TTLs(:,1) == session, :); % select only the required part
    [Srise, ~] = getTrialTTLs(tsessions_TTLs, TTLPath);
    Srise = double(Srise); % contains sample onset trial
    AudRise = [];
    SomRise = [];

    % Align multimodal sessions to both stimuli
    if strcmp(stimuli_parameters.Par.Rec, 'SxA')
        AudRise = Srise; % align multimodal trials to sound onset
        sa_idx = stimuli_parameters.Stm.Var25==4 & ...
            stimuli_parameters.Stm.SomAudSOA < 0;
        AudRise(sa_idx) = Srise(sa_idx) + round(abs(stimuli_parameters.Stm.SomAudSOA(sa_idx)) ./ Fs);
        alignspikes(spiketimes, size(spiketimes, 1), AudRise, stimuli_parameters, Fs);

        SomRise = Srise; % align multimodal trials to som onset
        sa_idx = stimuli_parameters.Stm.Var25==4 & ...
            stimuli_parameters.Stm.SomAudSOA > 0;
        SomRise(sa_idx) = Srise(sa_idx) + round(abs(stimuli_parameters.Stm.SomAudSOA(sa_idx)) ./ Fs);
        aligned_spikes = alignspikes(spiketimes, size(spiketimes, 1), SomRise, stimuli_parameters, Fs);
    else
        aligned_spikes = alignspikes(spiketimes, size(spiketimes, 1), Srise, stimuli_parameters, Fs);
    end

    % save spike alignment
    filename = sprintf('M%.2i_S%.2i_%s_AlignedSpikes', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
    save(fullfile(OutPath, filename), "aligned_spikes","Srise","AudRise","SomRise")
    clearvars("sa_idx")
end


%% align spikes
% saves aligned spikes with corresponding stimulus onset and offset timings
% spike times in sec, Srise & Sfall in samples

Copy_of_alignspikes(BehaviorPath, TTLPath, OutPath, spiketimes, relevant_sessions, skip_sessions, cids, sessions_TTLs, Fs);

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

%% optional: rename mouse nr

clearvars

BehaviorPath = 'D:\DATA\Behavioral Stimuli\M3\'; % stimuli parameters
SoundPath = [BehaviorPath 'Sound recording\'];
%OutPath = 'D:\DATA\Processed\M23\Spectograms';
% make array with all relevant session numbers
sound_file_list = dir(fullfile(SoundPath, '*_Sound.mat'));
stim_file_list = dir(fullfile(BehaviorPath, '*.mat'));

% load behaviour file
for session_file = 1:size(stim_file_list,1)
    %sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
    %stim_file = dir(fullfile(BehaviorPath, sessionFile));
    % load files
    stimuli_parameters = load([stim_file_list(session_file).folder '\' stim_file_list(session_file).name], 'Stm', 'Par');
    Stm = stimuli_parameters.Stm;
    Par = stimuli_parameters.Par;

    % rename in Par file
    Par.MouseNum = num2str(23);

    % rename in Stm file
    Stm.MouseNum(:) = 23;

    save([stim_file_list(session_file).folder '\' stim_file_list(session_file).name], 'Stm', 'Par'); % stimuli_parameters

end