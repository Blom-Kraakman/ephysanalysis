%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   plot single unit responses & single units in channel map

clearvars

% set directories
recordingFolder = 'D:\DATA\EphysRecordings\M9\M09_2024-04-23_12-46-03\Record Node 103\experiment1\recording1\';
recPath = [recordingFolder 'continuous\Intan-100.Rhythm Data-A\'];
TTLPath = [recordingFolder 'events\Intan-100.Rhythm Data-A\TTL\'];
messagesPath = [recordingFolder 'events\MessageCenter\'];
KSPath = 'D:\DATA\EphysRecordingsSorted\M09_2\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M9\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M9_2'; % output directory

rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording

%relevant_sessions = [1 9]; % M6 behaviour files (if only 1 behavior file in rec: [1 1])
relevant_sessions = [4 11]; % M9
skip_sessions = 0;
%relevant_sessions = [1 11]; % M8
%skip_sessions = 10; % M8

%relevant_sessions = [1 10]; %M11 and M10
%skip_sessions = 2;

Fs = 30000; % sampling freq


%% sessions TTLs
% save session TTLs if needed, else load correct file

TTLs_file = dir([OutPath '\*_OE_TTLs.mat']);

if isempty(TTLs_file)
    [sessions_TTLs, ~] = getSessionTTLs(messagesPath);
    filename = sprintf('M%.2i_S%02d-%02d_OE_TTLs', str2double(Par.MouseNum), relevant_sessions(1), relevant_sessions(2));
    save(fullfile(OutPath, filename), "sessions_TTLs")
    disp('sessions TTLs saved')
else
    TTLs = load([OutPath '\' TTLs_file.name]); 
    sessions_TTLs = TTLs.sessions_TTLs;
    disp('sessions TTLs loaded from saved file')
end

%% Kilosort: post-curation unit extraction
%IronClust: post-curation unit extraction [spiketimes, cids,cpos] = ircGoodClusters(spiketimecsv,clusterqualitycsv);
% spike extraction from curated units

[spiketimes, cids, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, TTLPath, relevant_sessions, skip_sessions, rec_samples, sessions_TTLs, Fs, OutPath);
%% remove double spikes from originating in Phy
% make into function
% get all spiketimes from each good unit
spiketimes2 = cell(length(cids), 1);

for cluster = 1:length(cids)
    Tspkt = spiketimes{cluster};
    ISI = zeros(length(Tspkt),1);

    %find minimal time between two spikes
    for spike = 1:(length(Tspkt)-1)
        ISI(spike) = min(abs(Tspkt(spike) - Tspkt(spike+1)));
    end

    % select spikes following next one with more than 0.2ms (6 samples)
    index = ISI>=6;
    spiketimes2{cluster} = Tspkt(index);
    disp(size(spiketimes{cluster}, 1) - size(spiketimes2{cluster}, 1))
end

spiketimes = spiketimes2;
%% align spikes
% alignment of extracted spikes to stimulus on/off-set
% spike times in sec

alignspikes(BehaviorPath, OutPath, spiketimes, relevant_sessions, skip_sessions, Srise, Sfall, cids, Fs);

%% FRA analysis
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster
close all
% select correct input files
aligned_spikes = load([OutPath, '\M09_S05_FRA_AlignedSpikes']);
stimuli_parameters = load([BehaviorPath 'M9_S05_FRA.mat']);

% function saves figures, change mouse name
FSL = 0;
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
session = 2;

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


%% plot single unit
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M10\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M10'; % output directory
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

cluster = 1; % specify cluster position
session = 2;

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
aligned_spikes = aligned_spikes.SpkT;

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

% add spacing where needed
idx = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
for ii = idx'
    aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} + 0.25;
end

% define shared x-lim parameters
preT  = -0.2;
postT = 1.5;
xrange = [preT, postT];
binsize = 0.01;
start_aud = 0;
start_som = 0.25;
end_som = 0.75;
end_aud = 1;

xlinerange = [start_aud start_som end_som end_aud];

fig = figure;
ax = gca;

index = stimuli_parameters.Stm.Amplitude == 0.3 & stimuli_parameters.Stm.Var25 == 2;

SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
StimType = stimuli_parameters.Stm.Var25(index);
Var = [StimType, SOM_Hz, SOM_Amp];
%Var =  [stimuli_parameters.Stm.Var25, stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude];
raster_yinc = [5,10,10];
raster_color = [0, 0, 0];

[f, YTick, ~, ~, ~, YTickLim] = plotraster(ax, aligned_spikes(:, cluster), Var, raster_color, raster_yinc, 1);

yticks(YTick{1});
yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
ylim(f,yrange);
xlim(f,xrange);
xline(xlinerange) % on/off set

for i = 1:size(YTickLim,1) % delimit groups
    yline(ax,YTickLim(i,1)-3,':k');
    yline(ax,YTickLim(i,2)+3,':k');
end

xlabel('Time (s)')
all_freqs = unique(stimuli_parameters.Stm.SomFreq);
yaxislabels = all_freqs(2:9);    
%yaxislabels = {'Control', 'Vibrotactile only', 'Vibrotactile & noise' 'Noise only'};

yticklabels(yaxislabels)


%% plotting SOM sessions
%saveplots = 0; %0 don't save, 1 save plots in OutPath
% currently plots raster + psth of each vibration freq seperatly

close all
OutPath = 'D:\DATA\Processed\M10'; % output directory

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

sessions = 2; % [exp ctrl]
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

fig = SOMplotting(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, 0);
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
plot_in_channel_map(KSPath, cpos);

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
