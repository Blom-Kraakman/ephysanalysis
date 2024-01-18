%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   quatify which cluster responds to which event
%       response = fire (spiketime) in event window
%   quantify firing rates? changes expected after events vs spontaneous

clearvars

% set directories
recPath = 'D:\DATA\EphysRecordings\M4\M04_2023-12-14_12-47-41\Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data\';
TTLPath = 'D:\DATA\EphysRecordings\M4\M04_2023-12-14_12-47-41\Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data\TTL\';
messagesPath = 'D:\DATA\EphysRecordings\M4\M04_2023-12-14_12-47-41\Record Node 103\experiment1\recording1\events\MessageCenter\'; % session TTLs
KSPath = 'D:\DATA\EphysRecordingsSorted\M04\trimmed\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M4\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M4'; % output directory

rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording

relevant_sessions = [1 23]; % behaviour files (if only 1 behavior file in rec: [1 1])

%% IronClust: post-curation unit extraction

%[spiketimes, cids,cpos] = ircGoodClusters(spiketimecsv,clusterqualitycsv);

%% Kilosort: post-curation unit extraction
% spike extraction from curated units

[spiketimes, cids, cpos, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, TTLPath, messagesPath, relevant_sessions, rec_samples);

% save details good units
set = sprintf('%02d-%02d', relevant_sessions(1), relevant_sessions(2));
filename = ['M04_S' set '_InfoGoodUnits'];
save(fullfile('D:\DATA\Processed', filename), "cpos") %cpos variables: unit id, channel, depth, avg firing rate, nr spikes;

%% align spikes
% alignment of extracted spikes to stimulus on/off-set
% spike times in sec
% good to have: TTL signaling start & end of stim presentation set

% function saves aligned spikes cell array, change mouse name
alignspikes(BehaviorPath, spiketimes, relevant_sessions, Srise, Sfall, cids);

%% FRA analysis
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster

% select correct input files
aligned_spikes = load('D:\DATA\Processed\M04_S01_FRA_AlignedSpikes');
stimuli_parameters = load([BehaviorPath 'M4_S01_FRA.mat']);

% function saves figures, change mouse name
FRAanalysis(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, 1);

%% plotting single sessions
% FRA: raster
% AMn: raster + PSTH
% SOM: raster + PSTH

%cids = load('D:\DATA\Processed\M04_S01-23_InfoGoodUnits.mat'); 

% select which session to plot
session = 22;

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile('D:\DATA\Processed', sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);

plotResponses(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath);

%% plotting SOM sessions with their controls
%saveplots = 0; %0 don't save, 1 save plots in OutPath
sessions = [23 21]; % [exp ctrl]
%cids = load('D:\DATA\Processed\M04_S01-23_InfoGoodUnits.mat'); 
cids = [    90   111   124   126   151   159   169   171   182   188   218   232   256   261   264   265   268   278];
OutPath = 'D:\DATA\Processed\M4';

SOMplotting(sessions, cids, OutPath, BehaviorPath, 0);
%% FSL SOM analysis
% function SOM = SOManalysis(stimuli_parameters, aligned_spikes, cids)
% input: stimuli_parameters.Par, stimuli_parameters.Stm, aligned_spikes
% output: first spike latency SOM/AM trials

% select which session to plot
session = 23;

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile('D:\DATA\Processed', sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
aligned_spikes = aligned_spikes.SpkT;

cids = [90 111 124 126 151 159 169 171 182 188 218 232 256 261 264 265 268 278];

SOM = FSL_SOM_AMn(stimuli_parameters, aligned_spikes, cids);

%% plotting channel map
% match unit position to channel map

channel_map = readNPY([KSPath 'channel_map.npy']);
channel_positions = readNPY([KSPath 'channel_positions.npy']);
cids = load('D:\DATA\Processed\M04_S01-23_InfoGoodUnits.mat'); 
cpos = cids.cpos;

xcoords = channel_positions(:,1);
ycoords = channel_positions(:,2);
margin = 50;
xMin = min(xcoords); xMax = max(xcoords);
xSpan = xMax-xMin; xMid = 0.5*(xMax+xMin);
yMin = min(ycoords); yMax = max(ycoords);
ySpan = yMax-yMin; yMid = 0.5*(yMax+yMin);
maxSpan = max(xSpan,ySpan);

fig = scatter(xcoords, ycoords, ".", 'k');
hold on;

for i = 1:length(channel_map)
    text((channel_positions(i,1)+1), (channel_positions(i,2)+1),num2str(channel_map(i)))
end

% add position of units included in analysis
hold on
units = table2array(cpos(:,2));
idx = ismember(channel_map, units);
fig = scatter(channel_positions(idx,1), channel_positions(idx,2), 'o'); % wrongly indexed

% to do: reflect number of clusters on 1 channel. apply jitter in circle
% marking to achieve this
% to do: add unit number with corresponding channel
%text((channel_positions(unique(idx),1)+1), (channel_positions(unique(idx),1)+1),num2str(channel_map()))

% format figure
title('Channel map')
ylabel('Relative depth')
fig.SizeData = 100;
fontsize(13,"points")
axis square
xlim([xMid - 0.5*maxSpan - margin, xMid + 0.5*maxSpan + margin]);
ylim([yMid - 0.5*maxSpan - margin, yMid + 0.5*maxSpan + margin]);
xticks(unique(xcoords));
yticks(unique(ycoords));


%% channel waveforms
% extract and plot waveform traces
addpath('C:\Users\TDT\Documents\GitHub\ephys_analysis');
for k=1:length(cids)
    close all
    %F = myfig(1,'fig');
    F = figure;
    F.Position = [2405,214,838,890];
    F = ChannelWaveforms(F,KSPath,cids(k),cpos(k,3),200);
    F.Renderer = 'painter';
    F.PaperUnits = 'inches';
    F.PaperSize = F.Position(3:4)./96;
    %saveas(F,[UserPath 'Processed\' DirName '\Units\' DirName '-U' num2str(cids(k)) '.pdf']);
    %     saveas(F,[UserPath 'Processed\' Mouse '\Units\' Mouse '-U' num2str(cids(k)) '.png']);
end

%% quantify reactive units
% to do

