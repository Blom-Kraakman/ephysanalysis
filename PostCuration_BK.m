%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   plot single unit responses & single units in channel map

clearvars

% set directories
recordingFolder = 'D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 103\experiment1\recording1\';
recPath = [recordingFolder 'continuous\Intan-100.Rhythm Data-A\'];
TTLPath = [recordingFolder 'events\Intan-100.Rhythm Data-A\TTL\'];
messagesPath = [recordingFolder 'events\MessageCenter\'];
KSPath = 'D:\DATA\EphysRecordingsSorted\M08\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M8\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M8'; % output directory

rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording

relevant_sessions = [1 11]; % M8
% relevant_sessions = [1 9]; % M6 behaviour files (if only 1 behavior file in rec: [1 1])
skip_sessions = 0;
Fs = 30000; % sampling freq

%% sessions TTLs as extracted from OpenEphys message center

[sessions_TTLs, sessions_TTLs_details] = getSessionTTLs(messagesPath, rec_samples, Fs);

%% save session TTLs
set = sprintf('%02d-%02d', relevant_sessions(1), relevant_sessions(2));
filename = ['M06_S' set '_OE_TTLs'];
save(fullfile(OutPath, filename), "sessions_TTLs")
%filename = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
filename = sprintf('M%.2i_S%02d-%02d_OE_TTLs', str2double(stimuli_parameters.Par.MouseNum), relevant_sessions(1), relevant_sessions(2));

%% Kilosort: post-curation unit extraction
%IronClust: post-curation unit extraction [spiketimes, cids,cpos] = ircGoodClusters(spiketimecsv,clusterqualitycsv);
% spike extraction from curated units

[spiketimes, cids, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, TTLPath, messagesPath, relevant_sessions, skip_sessions, rec_samples, Fs, OutPath);

%% align spikes
% alignment of extracted spikes to stimulus on/off-set
% spike times in sec

alignspikes(BehaviorPath, OutPath, spiketimes, relevant_sessions, skip_sessions, Srise, Sfall, cids, Fs);

%% FRA analysis
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster
close all
% select correct input files
aligned_spikes = load([OutPath, '\M08_S08_FRA_AlignedSpikes']);
stimuli_parameters = load([BehaviorPath 'M8_S08_FRA.mat']);

% function saves figures, change mouse name
FSL = 1;
FRAanalysis(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, FSL);

%% plotting single sessions
% FRA: raster
% AMn: raster + PSTH
% SOM: raster + PSTH

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

% select which session to plot
%close all
session = 4;

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);

if strcmp(stimuli_parameters.Par.Rec, 'SxA')
    idx = strcmp(stimuli_parameters.Stm.MMType, "SO");
    stimuli_parameters.Stm(idx,25) = {2};
    idx = strcmp(stimuli_parameters.Stm.MMType, "SA");
    stimuli_parameters.Stm(idx,25) = {3}; 
    idx = strcmp(stimuli_parameters.Stm.MMType, "OA");
    stimuli_parameters.Stm(idx,25) = {4}; 
    idx = strcmp(stimuli_parameters.Stm.MMType, "OO");
    stimuli_parameters.Stm(idx,25) = {1};
    % order: type, freq, amplitude
end

plotResponses(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath);

%% plotting SOM sessions
%saveplots = 0; %0 don't save, 1 save plots in OutPath
% currently plots raster + psth of each vibration freq seperatly

close all
OutPath = 'D:\DATA\Processed\M8'; % output directory

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

sessions = 4; % [exp ctrl]
%load correct files
session = sessions(1); % select experimental session

sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);

% plot with matching control
% session = sessions(2); % select control session
% 
% sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
% stim_files = dir(fullfile(BehaviorPath, sessionFile));
% stimuli_parameters_ctrl = load([stim_files.folder '\' stim_files.name]);
%
% aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
% aligned_spikes_ctrl = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
% aligned_spikes_ctrl = aligned_spikes_ctrl.SpkT;
% %format data to plot together
% stimuli_parameters = vertcat(stimuli_parameters_som.Stm, stimuli_parameters_ctrl.Stm);
% aligned_spikes = vertcat(aligned_spikes_som, aligned_spikes_ctrl);

OutPath = 'D:\DATA\Processed\M8\test'; % output directory

fig = SOMplotting(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, 1);
%% FSL SOM analysis
% function SOM = SOManalysis(stimuli_parameters, aligned_spikes, cids)
% input: stimuli_parameters.Par, stimuli_parameters.Stm, aligned_spikes
% output: first spike latency SOM/AM trials

% select which session to plot
session = 2;

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
aligned_spikes = aligned_spikes.SpkT;


%% maybe try swarmchart
SOM = FSL_SOM_AMn(stimuli_parameters, aligned_spikes, cids);

%% plotting channel map
% match unit position to channel map

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos;
plot_in_channel_map(KSPath, cids);

%% channel waveforms
% extract and plot waveform traces
%addpath('C:\Users\TDT\Documents\GitHub\ephys_analysis');
%cids = [90 111 124 126 151 159 169 171 182 188 218 232 256 261 264 265 268 278];
cids = load([OutPath, '\M04_S01-23_InfoGoodUnits.mat']); 
cpos = cids.cpos;

for k=1:size(cpos, 1)
    close all
    %F = myfig(1,'fig');
    F = figure;
    F.Position = [2405,214,838,890];
    F = ChannelWaveForms(F, KSPath, cpos.id(k), cpos.channel(k), 200);
    F.Renderer = 'painter';
    F.PaperUnits = 'inches';
    F.PaperSize = F.Position(3:4)./96;
    %saveas(F,[UserPath 'Processed\' DirName '\Units\' DirName '-U' num2str(cids(k)) '.pdf']);
    %     saveas(F,[UserPath 'Processed\' Mouse '\Units\' Mouse '-U' num2str(cids(k)) '.png']);
end

%% plot spikes of whole session
% to do

%% plot vibrotac stimulus signal
% include also feedback signal
% see plotting in GUI
% not correct yet, Ramp and offset

Waveform = 'BiSine';
Fs = 30000;
StimDur = 500;
Amplitude = 1;
SomFreq = 100;
Ramp = 10;
ISI = 1000;
Offset = 0.1;

StimDurSamp = ceil(StimDur * 0.001 * Fs);
som_waveform = nan(1,StimDurSamp);
tt = 0:1/Fs:(StimDur* 0.001);
switch Waveform
    case {'Square'}
        som_waveform = Amplitude .* ones(1,StimDurSamp);
        som_waveform(1) = 0; som_waveform(end) = 0; % zero at beginning or end
    case {'UniSine'}
        som_waveform =  0.5 * Amplitude .* ( 1-cos(2*pi*SomFreq .* tt) );
    case {'BiSine'}
        som_waveform =  Amplitude .* ( sin(2*pi*SomFreq .* tt) );
end

% apply envelope (On-/Off-ramps)
if Ramp > 0
    Nenv			=	round( Ramp * 10^-3 * Fs );
    som_waveform    =	envelope(som_waveform',Nenv)';
end

% padding pre-post stimulus time with zeroes
PrePostDur = 0.2 * ISI * 0.001;
PrePostSamp = round(PrePostDur * Fs);
tt = [-flip(1:PrePostSamp)./Fs, tt, tt(end) + (1:PrePostSamp)./Fs];
som_waveform = [zeros(1,PrePostSamp), som_waveform, zeros(1,PrePostSamp)];

% plotting the waveform
plot(tt(1:length(som_waveform)),som_waveform+Offset);
xlim([min(tt),max(tt)]);
