%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   plot single unit responses & single units in channel map

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
%% remove double spikes from originating in Phy - No longer needed
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

%spiketimes = spiketimes2;

%% OLD get trials onset TTLs of all sessions in recording OLD
[Srise, Sfall] = TTLsToUse(sessions_TTLs, TTLPath, rec_samples);

%% align spikes
% alignment of extracted spikes to stimulus on/off-set
% spike times in sec

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

%% plotting single sessions
% FRA: raster
% AMn: raster + PSTH
% SOM: raster + PSTH
close all

OutPath = 'D:\DATA\Processed\M14'; % output directory
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M14\'; % stimuli parameters

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

% select which session to plot
%session = 18;
for session = 1:5
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
end

%% plot single unit
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M15\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M15\ICC'; % output directory
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

cluster = 1; % specify cluster position M10:441=20  M8 331=17
session = 17; % m10 2; m8 5

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
aligned_spikes = aligned_spikes.SpkT;

if strcmp(stimuli_parameters.Par.Rec, 'SxA')
    idx = strcmp(stimuli_parameters.Stm.MMType, "SO");
    stimuli_parameters.Stm(idx,25) = {3};
    idx = strcmp(stimuli_parameters.Stm.MMType, "SA");
    stimuli_parameters.Stm(idx,25) = {2}; 
    idx = strcmp(stimuli_parameters.Stm.MMType, "OA");
    stimuli_parameters.Stm(idx,25) = {4}; 
    idx = strcmp(stimuli_parameters.Stm.MMType, "OO");
    stimuli_parameters.Stm(idx,25) = {1};
    % order: type, freq, amplitude
end

% % add spacing where needed
% idx = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
% for ii = idx'
%     aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} + 0.25;
% end

% define shared x-lim parameters
preT  = -0.2;
postT = 1.2;
xrange = [preT, postT];
binsize = 0.01;
start_aud = 0;
start_som = 0.25;
end_som = 0.75;
end_aud = 1;

xlinerange = [start_aud start_som end_som end_aud];

%fig = figure; % rasterplot
fig = subplot(2,1,1); % rasterplot
ax = gca;

index = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600);

SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
StimType = stimuli_parameters.Stm.Var25(index);
Var = [SOM_Hz, StimType];
%Var =  [stimuli_parameters.Stm.Var25, stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude];
raster_yinc = [5,10,10];
raster_color = [0, 0, 0];

[f, YTick, ~, ~, ~, YTickLim] = plotraster(ax, aligned_spikes(index, cluster), Var, raster_color, raster_yinc, 1);

yticks(YTick{1});
yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
ylim(f,yrange);
xlim(f,xrange);
xline(xlinerange) % on/off set

for i = 1:size(YTickLim,1) % delimit groups
    yline(ax,YTickLim(i,1)-3,':k');
    yline(ax,YTickLim(i,2)+3,':k');
end

%xlabel('Time (s)')
all_freqs = unique(stimuli_parameters.Stm.SomFreq);
yaxislabels = {'Control',all_freqs(2:8)};
yticklabels(yaxislabels)
ylabel('Vibrotactile frequency (Hz)')
ax.FontSize = 16;

% make histogram / PSTH
binsize = 0.02;

fig = subplot(2,1,2); 
hold on
% select groups for hist
OO = stimuli_parameters.Stm.Var25 == 1;
OA = stimuli_parameters.Stm.Var25 == 4;
SO = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600) & stimuli_parameters.Stm.Var25 == 3;
SA = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600) & stimuli_parameters.Stm.Var25 == 2;

[N, edges] = histcounts(vertcat(aligned_spikes{OO, cluster}), preT:binsize:postT); % OO
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(OO))/binsize)) % spike/s
%plot(edges(1:end-1), ((N/sum(OO))/binsize), 'Color', '#D95319', 'LineWidth',1.5) % spike/s
plot(edges(1:end-1), ((N/sum(OO))/binsize), 'k','LineWidth',2) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{OA, cluster}), preT:binsize:postT); % OA
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(OA))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(OA))/binsize), 'Color', "#CE87AA",'LineWidth',2) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{SO, cluster}), preT:binsize:postT); % SO
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(SO))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(SO))/binsize), 'Color', "#8BC9E8",'LineWidth',2) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{SA, cluster}), preT:binsize:postT); % SA
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(SA))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(SA))/binsize), 'Color', "#9474A6",'LineWidth',2) % spike/s

%format axis
xlabel('Time (s)')
ylabel('Spike rate (spikes/s)')
xline(xlinerange) % on/off set
legend('control', 'sound only', 'vibrotactile only', 'multimodal', 'Location', 'northeast')

ax2 = gca;
ax2.FontSize = 16;
xlim(ax2,xrange);
ylim(ax2, [0 100])


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

fig = SOMplotting(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, 0);

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
