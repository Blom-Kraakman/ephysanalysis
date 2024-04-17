%% Quantify results
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   quatify which cluster responds to which event
%       response = fire (spiketime) in event window
%   quantify firing rates? changes expected after events vs spontaneous

clearvars

% set directories
recordingFolder = 'D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\';
%acfeedbackPath = [recordingFolder 'Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B'];
recPath = [recordingFolder 'Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\'];
TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
%messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
%rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
KSPath = 'D:\DATA\EphysRecordingsSorted\M08\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M8\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M8'; % output directory

Fs = 30000; % sampling freq

% select which session to analyse
session = 4;

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
if isfield(aligned_spikes,"Srise")
    disp("Srise/Sfall loaded from data file.")
    Srise = aligned_spikes.Srise;
    Sfall = aligned_spikes.Sfall;
end
aligned_spikes = aligned_spikes.SpkT;

% load sessions details
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']).name;
sessions_TTLs = load([OutPath '\' TTLs_file]);

% Kilosort: post-curation unit extraction
% needed to get Srise...
rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
relevant_sessions = [1 11]; % M8
skip_sessions = 10; % M8
[spiketimes, cids, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, TTLPath, relevant_sessions, skip_sessions, rec_samples, sessions_TTLs.sessions_TTLs, Fs, OutPath);

%% select TTLs of session 4
nTrials = size(aligned_spikes, 1);
%nClusters = size(aligned_spikes, 2);
nClusters = 1;
NStim = 1;

% get TTL on and off for a session
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));
for file = 1:9

    stimuli_parameters_temp = load([stim_files(file).folder '\' stim_files(file).name]);

    disp(NStim);

    % session to plot
    if ismember(str2double(stimuli_parameters_temp.Par.Set), session)
        % keep Srise and Sfall withing boundaries
        tSrise = Srise(NStim: (NStim + size(stimuli_parameters_temp.Stm, 1)-1));
        tSfall = Sfall(NStim: (NStim + size(stimuli_parameters_temp.Stm, 1)-1));
    end

    % update cummulative stimuli
    NStim = NStim + size(stimuli_parameters_temp.Stm, 1);

end
clearvars("stimuli_parameters_temp")

%% get and plot analog signal

% read and save into variable
fid = fopen('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\continuous.dat');
gain = 10 / 2^16; % V/bit; estimated +/- 5V / 16-bit
analog_trace = gain*double(fread(fid, '*int16'));
analog_samples = readNPY('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\sample_numbers.npy');
% --- warning: assuming there are no pauses within the analog recording ---
tSriseIdx = tSrise - analog_samples(1) + 1;
tSfallIdx = tSfall - analog_samples(1) + 1;
%% ABW --- Start

% ~~~ select subset of data for testing ~~~
all_freqs = unique(stimuli_parameters.Stm.SomFreq);
amp = 0.1;
nClusters = length(cids);
cluster = 1; % cluster index
cluster = find(cids == 331); %find(cids == 331) % find(cids == 184)
freq = find(all_freqs == 10);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nCol = 4;
nRow = 4;

% Part 1: raster plot with analogue signal on top

% comment next lines for testing; uncomment for loop
% dont forgot the `end` statement at the end
% for cluster = 1:nClusters
    for freq = 1:length(all_freqs)
        figure;

        index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == amp);
        SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
        Var = SOM_Hz;
        yaxistext = [num2str(all_freqs(freq)) ' Hz'];

        % make rasterplot
        f = subplot(nRow,nCol,[nCol+1,nRow*nCol-1]);
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(f, aligned_spikes(index, cluster), Var, [0, 0, 0], [5 5], 1);
        yrange = [min(YTick{end}) - 50, max(YTick{end}) + 50];

        % get analogue signal
        PreT = 0.25; PostT = .75; %s
        PreT_samp = -round(PreT*Fs);
        PostT_samp = round(PostT*Fs);
        trialon  = double(tSriseIdx(index));
        trialoff = double(tSfallIdx(index));
        
        tt = (PreT_samp:PostT_samp) ./ Fs;
        motorSignal = double(analog_trace(trialon + (PreT_samp:PostT_samp)));
        ax2 = subplot(nRow,nCol,[1,nCol-1]);
        hold(ax2,"off");
        plot(ax2,tt,motorSignal','Color',[.5,.5,.5]);hold(ax2,"on");
        plot(ax2,tt,mean(motorSignal,1),'Color','k');hold(ax2,"off");
        linkaxes([f,ax2],'x');
        xlim(ax2,[-PreT,PostT])
        xlabel(f,"Time (s)")

% Part 2 cycle

% index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == amp);
% SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
Mf = unique(SOM_Hz);
% ~~~ parameters for analysis ~~
StimDur = 0.5; %s
% StimDur = 1/Mf; % 1/Mf = 1 cycle
SkipVal = 0;
SkipMethod = '';
nBins = 24;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% calculation
[CycT] = CycTimes(aligned_spikes(index, cluster),StimDur, Mf,SkipVal, SkipMethod);
[Vs,Ph,Ray] = calcVS(CycT);

% plotting
ax_cycHistPol = subplot(nRow,nCol,[nCol,nCol+nRow],polaraxes);
polarhistogram(ax_cycHistPol,2*pi*cat(1,CycT{:}),nBins);
hold(ax_cycHistPol,"on");
RMax = ax_cycHistPol.RLim(2);
plot(ax_cycHistPol,[0,2*pi*Ph],[0,Vs*RMax],'.r-','LineWidth',2)
ax_cycHistCat = subplot(nRow,nCol,[nCol+2*nRow,nCol+3*nRow]);
histogram(ax_cycHistCat,cat(1,CycT{:}),nBins)
xline(Ph,'r--','LineWidth',2)
ylabel(ax_cycHistCat,'Number of spikes')
xlabel(ax_cycHistCat,'Phase (cycle)')
xlim(ax_cycHistCat,[0,1])
title(ax_cycHistCat,...
    {['Vector strength: ',num2str(Vs,'%.2f')],...
    ['Phase: ',num2str(Ph,'%.2f cycle')],...
    ['p-value: ', num2str(Ray)]})



sgtitle([num2str(Mf),' Hz   ', num2str(amp),' V',' - ', 'unit ',num2str(cids(cluster))])
    %end
end

% ABW --- End