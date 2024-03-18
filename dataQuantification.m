%% Quantify results
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   quatify which cluster responds to which event
%       response = fire (spiketime) in event window
%   quantify firing rates? changes expected after events vs spontaneous

clearvars

% set directories
recordingFolder = 'D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\';
acfeedbackPath = [recordingFolder 'Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B'];
recPath = [recordingFolder 'Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\'];
%TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
%messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
%KSPath = 'D:\DATA\EphysRecordingsSorted\M08\'; % kilosort ephys data
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

% load sessions details
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']).name;
sessions_TTLs = load([OutPath '\' TTLs_file]);



%% calculate firing rates pre and post simulus onset

% onset stimuli
PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/4)/1000; % baseline period
PostT = 0.05; % post stim period

% calculate baseline firing rate
[baselineRate, stimulusRate] = baselinerate(aligned_spikes.SpkT, PreT, PostT); % get firing rates (Hz) in time window

%% plot summary fig - vibrotactile sessions
% boxplot, per unit

close all
amp = 0.3; 

uFreq = unique(stimuli_parameters.Stm.SomFreq);
nFreq = length(uFreq);
nClusters = length(cids);
condition = nan(nClusters, nFreq);
figure

for cluster = 1:nClusters
    for freq = 1:nFreq
        condition(cluster, freq) = median(stimulusRate([stimuli_parameters.Stm.SomFreq] == uFreq(freq) & [stimuli_parameters.Stm.Amplitude] == amp, cluster))';
    end

    plot(uFreq, mean(condition(cluster, :) - condition(cluster, 1), 1), ':o')
    hold on

end

xticks(uFreq)
xticklabels(uFreq(1:nFreq));
xlabel(['Vibrotactile stimulus (Hz), ' num2str(amp) ' (V)'])
ylabel('Normalized firing rate (Hz)')

title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

save figure
figname = sprintf('M%.2i_S%.2i_%s_mean', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

% plot summary fig: boxplot, all units
close all
condition = nan(nClusters, nFreq);
figure

% select data
for freq = 1:nFreq
    condition(:, freq) = median(stimulusRate([stimuli_parameters.Stm.SomFreq] == uFreq(freq) & [stimuli_parameters.Stm.Amplitude] == amp, :))';
end

% plot
figure
boxplot(condition)

xticklabels(uFreq(1:nFreq));
xlabel(['Vibrotactile stimulus (Hz), ' num2str(amp) ' (V)'])
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

% save figure
figname = sprintf('M%.2i_S%.2i_%s_boxplot_all units', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));


%% plot summary fig - pressure sessions
% boxplot, per  units
close all

% parameters
uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
nClusters = length(cids);
condition = nan(nClusters, nAmp);

figure
for cluster = 1:nClusters

    % select data
    for amplitude = 1:nAmp
        condition(cluster, amplitude) = median(stimulusRate([stimuli_parameters.Stm.Amplitude] == uAmp(amplitude), cluster))';
    end

    % plot
    plot(uAmp, mean(condition(cluster, :) - condition(cluster, 1), 1), ':o');
    hold on

end

xticks(uAmp)
xticklabels(uAmp(1:nAmp));
xlabel('Pressure (V)')
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

%save figure
figname = sprintf('M%.2i_S%.2i_%s_mean', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

% plot summary fig: boxplot, all units
close all
condition = nan(nClusters, nAmp);

% select data
for amplitude = 1:nAmp
    condition(:, amplitude) = median(stimulusRate([stimuli_parameters.Stm.Amplitude] == uAmp(amplitude), :))';
end

% plot
figure
boxplot(condition)

xticklabels(uAmp(1:nAmp));
xlabel('Pressure (V)')
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

%save figure
figname = sprintf('M%.2i_S%.2i_%s_boxplot_all units', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

   
%% quantify reactive units
% 1. quantify responsive units: sig diff firing rate between stim period vs baseline

% onset stimuli
PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/5)/1000; % baseline period
PostT = 0.05; % post stim period

uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
uFreq = unique(stimuli_parameters.Stm.SomFreq);
nFreq = length(uFreq);
nClusters = length(cids);

results.cids = cids;
results.conditions = uAmp;
results.pvalue = nan(nAmp, nClusters);
results.signrank = nan(nAmp, nClusters);
results.hvalue = nan(nAmp, nClusters);
results.responsive = [];

% stats test difference baseline vs stimulus period
for cluster = 1:nClusters

    % calculate baseline firing rate
    [baselineRate, stimulusRate] = baselinerate(aligned_spikes.SpkT, PreT, PostT); % get firing rates (Hz) in time window
   
    %signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    control = stimuli_parameters.Stm.Amplitude == 0;

    for condition = 1:nAmp
        index = stimuli_parameters.Stm.Amplitude == uAmp(condition);
        [p,h,stats] = signrank(stimulusRate(control,cluster), stimulusRate(index,cluster), 'alpha', 0.01); % different from control trials?
        results.pvalue(condition, cluster) = p;
        results.signrank(condition, cluster) = stats.signedrank;
        results.hvalue(condition, cluster) = h;
    end

    % unit responsive if sig diff for at least one condition
    if max(results.hvalue(:, cluster)) == 1
        results.responsive = [results.responsive, results.cids(cluster)];
    end

end


%% quantify reactive units
% 2. cross correlating single trials (KDF as in previous script), corr
% coeficient over control value to be sig

%% quantify reactive units
% 3. firing in relation to phase stimulus (try cycle histogram)
% get recorded command per session

% read and save into variable
fid = fopen('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\continuous.dat');
analog_trace = fread(fid, '*int16');
analog_samples = readNPY('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\sample_numbers.npy');

%% sessions_TTLs.sessions_TTLs(9, :) = [];

% select which session to analyse
session = 4;

% get start + end of session
idx = sessions_TTLs.sessions_TTLs(:,1) == session;
session_time = sessions_TTLs.sessions_TTLs(idx, 3);

% get first and last sample from session
start_session = analog_samples(analog_samples(:, 1) == session_time(1));
end_session = analog_samples(analog_samples(:, 1) == session_time(2));

plot(analog_samples(start_session:end_session), analog_trace(start_session:end_session)); % analog trace during session 4

% define start/end cycle + cycle window
trialDur = ((str2double(stimuli_parameters.Par.SomatosensoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI))/1000) * Fs; %in samples
nr_cycles = (stimuli_parameters.Stm.SomFreq .* (str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000));
cycle_window = ((str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000) ./ nr_cycles); %samples

% get samples corresponding to each cycle window
% for each cell
for cluster = 1:size(aligned_spikes.SpkT, 2)
    % for all trials
    for stim = 1:size(aligned_spikes.SpkT, 1)
        % time window each cycle
        cycle_on = 0;
        cycle_off = cycle_window(stim);
        for cycle = 1:nr_cycles(stim)
            % get spike times in this time window from aligned_spikes
            taligned_spikes = aligned_spikes.SpkT{stim, cluster};
            idx = (taligned_spikes >= cycle_on) & (taligned_spikes <= cycle_off);
            tcycle_aligned{(stim)} = taligned_spikes(idx);
            
            % update window
            cycle_on = cycle_on + cycle_window(stim);
            cycle_off = cycle_off + cycle_window(stim);
        end

        cycle_aligned = [cycle_aligned, tcycle_aligned];
    end

end


% match matrix with firing rate



%% calculate baseline firing rate
function [baseFR, expFR] = baselinerate(SpkT, PreT, PostT)

nClust = size(SpkT,2);
nTrials = size(SpkT,1);
baseFR = nan(nTrials, nClust);
expFR = nan(nTrials, nClust);

for cluster = 1:size(SpkT,2)
    for trial=1:size(SpkT,1)
        
        % get spike times per trial
        S = SpkT{trial, cluster};
        
        % count spikes in time windows
        baseSpikes = sum(S <= 0);
        expSpikes = sum(S > 0 & S <= PostT);

        % convert to firing rate (Hz)
        baseFR(trial, cluster) = baseSpikes/abs(PreT);
        expFR(trial, cluster) = expSpikes/abs(PostT);
    end
end
end
