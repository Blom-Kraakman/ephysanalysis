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
%TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
%messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
%rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
%KSPath = 'D:\DATA\EphysRecordingsSorted\M08\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M8\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M8'; % output directory

Fs = 30000; % sampling freq

%% get data from session to analyse
session = 5;

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

% define stimulus onset, pre and post simulus onset

% onset stimuli
PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/4)/1000; % baseline period
PostT = 0.25; % post stim period % 0.25 for sxa

if strcmp(stimuli_parameters.Par.Rec, 'SxA')
            % add delay to "SO"
            StimOn = stimuli_parameters.Stm.SomAudSOA./1000;
            StimOn(isnan(stimuli_parameters.Stm.SomAudSOA)) = 0;
else
        StimOn = zeros(size(stimuli_parameters.Stm, 1), 1);
end

% calculate firing rate (Hz) in time window
[baselineRate, stimulusRate] = firingrate(aligned_spikes, PreT, PostT, StimOn);
   
%% quantify reactive units
% 1. quantify responsive units: sig diff firing rate between stim period vs baseline

if strcmp(stimuli_parameters.Par.Rec, 'AMn')
    uFreq = unique(stimuli_parameters.Stm.Mf);
    uAmp = unique(stimuli_parameters.Stm.Intensity);
else
    uAmp = unique(stimuli_parameters.Stm.Amplitude);
    uFreq = unique(stimuli_parameters.Stm.SomFreq);
end

nAmp = length(uAmp);
nFreq = length(uFreq);

nClusters = length(cids);

results.session = stimuli_parameters.Par.Set;
results.type = stimuli_parameters.Par.Rec;
results.cids = cids;
results.ampitudes = uAmp;
results.frequencies = uFreq;
results.pvalue = nan(nAmp, nFreq, nClusters);
results.signrank = nan(nAmp, nFreq, nClusters);
results.hvalue = nan(nAmp, nFreq, nClusters);
results.responsive = [];

% stats test difference baseline vs stimulus period
for cluster = 1:nClusters

    signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    for freq = 1:nFreq
        for condition = 1:nAmp

            % session specific parameters
            if strcmp(stimuli_parameters.Par.Rec, 'AMn') % noise & AM
                index = stimuli_parameters.Stm.Mf == uFreq(condition);
                baseline = baselineRate(index,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') %vibrotac
                control = (stimuli_parameters.Stm.Amplitude == uAmp(condition)) & (stimuli_parameters.Stm.SomFreq == 0);
                index = (stimuli_parameters.Stm.Amplitude == uAmp(condition)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') % pressure
                    control = stimuli_parameters.Stm.Amplitude == 0;
                    index = stimuli_parameters.Stm.Amplitude == uAmp(condition);
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                end
            elseif strcmp(stimuli_parameters.Par.Rec, 'SxA') % multimodal
                control = strcmp(stimuli_parameters.Stm.MMType, 'OO');
                %control = (stimuli_parameters.Stm.Amplitude == 0) & (stimuli_parameters.Stm.SomFreq == 0);
                index = (strcmp(stimuli_parameters.Stm.MMType, 'SO')) & (stimuli_parameters.Stm.Amplitude == uAmp(condition)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            end

            %[p,h,stats] = signrank(stimulusRate(control,cluster), stimulusRate(index,cluster), 'alpha', 0.01); % different from control trials?
            % stat test
            [p,h,stats] = signrank(baseline, stim_evoked, 'alpha', 0.01); % different from control trials?
            results.pvalue(condition, freq, cluster) = p;
            results.signrank(condition, freq, cluster) = stats.signedrank;
            results.hvalue(condition, freq, cluster) = h;
        end
    end

    % unit responsive if sig diff for at least one condition
    if (strcmp(stimuli_parameters.Par.Rec, 'SxA')) && (max(max(results.hvalue(2:nAmp, 2:nFreq, cluster))) == 1) % && (results.hvalue(1, 1, cluster) == 0)
        results.responsive = [results.responsive, results.cids(cluster)];
    elseif (strcmp(stimuli_parameters.Par.Rec, 'AMn')) && max(results.hvalue(:, :, cluster)) % && (results.hvalue(1, 1, cluster) == 0))
        results.responsive = [results.responsive, results.cids(cluster)];
    elseif (strcmp(stimuli_parameters.Par.Rec, 'SOM'))
        if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') && max(results.hvalue(:, :, cluster))
            results.responsive = [results.responsive, results.cids(cluster)];
        elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') && max(max(results.hvalue(:, :, cluster)))
            results.responsive = [results.responsive, results.cids(cluster)];
        end
    end
end

%responsive units contains all
responsive_units(str2double(stimuli_parameters.Par.Set)) = results;
resp_cids = find(ismember(responsive_units(session).cids, responsive_units(session).responsive)); % position of responsive units in cids

%% quantify response strength per condition

conditions = {"OO", "OA", "SO", "SA"};

for i = 1:length(conditions)
    index = strcmp(stimuli_parameters.Stm.MMType, conditions{i});
    firing_mean(i, :) = mean(stimulusRate(index, resp_cids));
    firing_sd(i,:) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
end

%plot bar graph
data = mean(firing_mean, 2);
errors = mean(firing_sd, 2);

figure;
bar(mean(firing_mean, 2)) % avg FR of all responsive units
hold on
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
%yline(data(4))
%scatter(1:length(data), firing_mean, 10, 'k', 'filled');
xlabel('conditions')
xticklabels(["OO", "OA", "SO", "SA"])
ylabel('mean stimulus evoked firing rate (Hz)')
title(PostT)

hold off

% normalized to control trials
data = data - data(1);
figure;
bar(data) % avg FR of all responsive units
xlabel('conditions')
xticklabels(["OO", "OA", "SO", "SA"])
ylabel('mean stimulus evoked firing rate (Hz)')

%% compare firing rate unimodal to multimodal stimulus presentation
% makes scatter plot of firing rates per responsive unit between two conditions

% calculate stimulus induced change in firing rate
dstimulusRate = stimulusRate - baselineRate;

conditions = {"OO", "OA", "SO", "SA"};
dfiring = nan(length(conditions), size(dstimulusRate, 2));

% for i = 1:length(conditions)
%     index = strcmp(stimuli_parameters.Stm.MMType, conditions{i});
%     dfiring(i, :) = mean(dstimulusRate(index, :));
% end

for i = 1:length(conditions)
    index = strcmp(stimuli_parameters.Stm.MMType, conditions{i});
    dfiring_mean(i, :) = mean(dstimulusRate(index, resp_cids));
    dfiring_sd(i,:) = std(dstimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
end

%plot
figure;
scatter(dfiring_mean(2, :), dfiring_mean(3,:)) % 1: SA, 2: SO
xlabel(['Firing changes during' conditions(2)])
ylabel(['Firing changes during' conditions(3)])
title(['Post stim window: ' num2str(PostT) 'ms'])

% which/how many cells prefer what?
somPref = [];
audPref = [];
noPref = [];

index = strcmp(stimuli_parameters.Stm.MMType, "OA");
dfiring_OA = dstimulusRate(index, resp_cids);
index = strcmp(stimuli_parameters.Stm.MMType, "SO"); % to do: select best freq
dfiring_SO = dstimulusRate(index, resp_cids);

for unit = 1:length(responsive_units(session).responsive) %length(dfiring_mean)

    %[p,h,stats] = signrank(stimulusRate(control,cluster), stimulusRate(index,cluster), 'alpha', 0.01); % different from control trials?
    % stat test
    [p,h,stats] = signrank(dfiring_OA(:, unit), dfiring_SO(:, unit), 'alpha', 0.01); % different from control trials?

    if h
        if abs(dfiring_mean(2,unit)) > abs(dfiring_mean(3,unit))
            audPref = [audPref; responsive_units(session).responsive(unit)];
        elseif abs(dfiring_mean(2,unit)) < abs(dfiring_mean(3,unit))
            somPref = [somPref; responsive_units(session).responsive(unit)];
        elseif dfiring_mean(2,unit) == dfiring_mean(3,unit)
            noPref = [noPref; responsive_units(session).responsive(unit)];
        end
    else
        noPref = [noPref; responsive_units(session).responsive(unit)];
    end

end


%%
figure;

bar([length(audPref) length(somPref)])

%% pressure tuning curve

% calculate stimulus induced change in firing rate
dstimulusRate = stimulusRate - baselineRate;

uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);

nClusters = length(cids);

figure;
hold on

for cluster = 1:nClusters

    for amp = 1:nAmp
        index = stimuli_parameters.Stm.Amplitude == uAmp(amp);
        toPlot(amp, cluster) = mean(dstimulusRate(index, cluster));
    end
    toPlot(2:5, cluster) = toPlot(2:5, cluster) - toPlot(1, cluster);

    if any(ismember(resp_cids, cluster))
        plot(toPlot(:, cluster), 'r');
        %resp_cells = [resp_cells; mean(toPlot(:, cluster))];
    else
        plot(toPlot(:, cluster), 'k');
    end

    ylabel('df Firing rate (Hz)')
    xlabel('stimulus strength (V)')
    xticks(1:nAmp)
    xticklabels(uAmp)
end

%% pressure tuning curve mean of responsive units
% calculate stimulus induced change in firing rate
dstimulusRate = stimulusRate - baselineRate;

uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);

nClusters = length(cids);

figure;
hold on

for amp = 1:nAmp
    index = stimuli_parameters.Stm.Amplitude == uAmp(amp);
    dfiring_mean(amp, :) = mean(stimulusRate(index, resp_cids));
    dfiring_sd(amp,:) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
end

plot(dfiring_mean, 'r')
plot(mean(dfiring_mean, 2), 'k', 'LineWidth', 1.25)
xticks(1:nAmp)
xticklabels(uAmp)

ylabel('df Firing rate (Hz)')
xlabel('stimulus strength (V)')

%% quantify reactive units
% 2. cross correlating single trials (KDF as in previous script), corr
% coeficient over control value to be sig

%% quantify reactive units
% 3. firing in relation to phase stimulus (try cycle histogram)
% get and plot analog signal

% read and save into variable
fid = fopen('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\continuous.dat');
analog_trace = fread(fid, '*int16');
analog_samples = readNPY('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B\sample_numbers.npy');

%% select which session to analyse
session = 4;

% get start + end of session
idx = sessions_TTLs.sessions_TTLs(:,1) == session;
session_time = sessions_TTLs.sessions_TTLs(idx, 3);

% get first and last sample from session
session_start = analog_samples(analog_samples(:, 1) == session_time(1));
session_end = analog_samples(analog_samples(:, 1) == session_time(end));

figure;
plot(analog_samples(session_start:session_end), analog_trace(session_start:session_end)); % analog trace during session 4

%% select TTLs of session 4
nTrials = size(aligned_spikes, 1);
%nClusters = size(aligned_spikes, 2);
nClusters = 1;
NStim = 1;

% get TTL on and off for a session
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));
for file = 1:9

    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);

    disp(NStim);

    % session to plot
    if ismember(str2double(stimuli_parameters.Par.Set), session)
        % keep Srise and Sfall withing boundaries
        tSrise = Srise(NStim: (NStim + size(stimuli_parameters.Stm, 1)-1));
        tSfall = Sfall(NStim: (NStim + size(stimuli_parameters.Stm, 1)-1));
    end

    % update cummulative stimuli
    NStim = NStim + size(stimuli_parameters.Stm, 1);

end
%% ABW --- Start
all_freqs = unique(stimuli_parameters.Stm.SomFreq);


% for cluster = 1:nClusters
cluster = 1;
    %for freq = 1:length(all_freqs)
freq = 3;
        figure;

        index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.1);
        SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
        Var = SOM_Hz;
        yaxistext = [num2str(all_freqs(freq)) ' Hz'];

        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(gca, aligned_spikes(index, cluster), Var, [0, 0, 0], [5 5], 1);
        yrange = [min(YTick{end}) - 50, max(YTick{end}) + 50];

        % get analogue signal
        PreT = 0.2; PostT = 1.5; %s
        PreT_samp = -round(PreT*Fs);
        PostT_samp = round(PostT*Fs);
        ttSrise = double(tSrise(index));
        ttSfall = tSfall(index);
        [trialon,~] = find(analog_samples == ttSrise');
        [trialoff,~] = find(analog_samples == ttSfall');

        tt = (PreT_samp:PostT_samp) ./ Fs;
        motorSignal = double(analog_trace(trialon + (PreT_samp:PostT_samp)));
        hold(f,'on');
        plot(f,tt,motorSignal','k');
    %end
% end

% [CycT] = CycTimes(aligned_spikes,StimDur, Mf,SkipVal, SkipMethod);


% ABW --- End
%% plot aligned spikes raster of each condition
all_freqs = unique(stimuli_parameters.Stm.SomFreq);

for cluster = 1:nClusters
    %for freq = 1:length(all_freqs)
freq = 7;
        figure;

        index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.1);
        SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
        Var = SOM_Hz;
        yaxistext = [num2str(all_freqs(freq)) ' Hz'];

        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(gca, aligned_spikes(index, cluster), Var, [0, 0, 0], [5 5], 1);
        yrange = [min(YTick{end}) - 50, max(YTick{end}) + 50];

    %end
end

%% plot analog motor signal per condition
for freq = 1:length(all_freqs)
index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.3);

% figure;
ttSrise = tSrise(index);
ttSfall = tSfall(index);

figure;
hold on

for i = 1:length(ttSrise)
    trialon = find(analog_samples == ttSrise(i));
    trialoff = find(analog_samples == ttSfall(i));
    tanalog_samples = trialon:trialoff;

    plot(double(analog_samples(tanalog_samples)-analog_samples(trialon))./Fs, analog_trace(tanalog_samples));

end

xlabel('Samples (normalized)')
title([num2str(all_freqs(freq)) 'Hz (session '  stimuli_parameters.Par.Set ')']);

end
%% define start/end cycle + cycle window
% does not work well

nCycles = (stimuli_parameters.Stm.SomFreq .* (str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000));
cycle_window = ((str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000) ./ nCycles); %samples
trialDur = ((str2double(stimuli_parameters.Par.SomatosensoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI))/1000) * Fs; %in samples

%cycle_aligned = nan(nClusters, nTrials, max(nCycles)); 
cycle_aligned = [];
% get samples corresponding to each cycle window
% for each cell
for cluster = 1:nClusters
    % for all trials
    for stim = 1:nTrials
        % time window each cycle
        cycle_on = 0;
        cycle_off = cycle_window(stim);
        for cycle = 1:nCycles(stim)
            % get spike times in this time window from aligned_spikes
            taligned_spikes = aligned_spikes.SpkT{stim, cluster};
            idx = (taligned_spikes >= cycle_on) & (taligned_spikes <= cycle_off);

            if max(idx) == 0; continue; end

            tcycle_aligned{(cycle)} = taligned_spikes(idx); % store spike times in cell array

            %cycle_aligned(cluster, stim, cycle) = taligned_spikes(idx);

            % update window
            cycle_on = cycle_on + cycle_window(stim);
            cycle_off = cycle_off + cycle_window(stim);
        end

        cycle_aligned = [cycle_aligned, tcycle_aligned];
    end

end


% match matrix with firing rate

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

%% calculate baseline firing rate
function [baseFR, expFR] = firingrate(SpkT, PreT, PostT, StimOn)

% SpkT: cell array with spike times aligned to stimulus onset
nClust = size(SpkT,2);
nTrials = size(SpkT,1);
baseFR = nan(nTrials, nClust);
expFR = nan(nTrials, nClust);

for cluster = 1:size(SpkT,2)
    for trial=1:size(SpkT,1)
        
        % get spike times per trial
        S = SpkT{trial, cluster};
        
        % count spikes in time windows
        baseSpikes = sum(S >= -PreT & S <= 0); % preT - 0
        expSpikes = sum(S > StimOn(trial) & S <= (StimOn(trial) + PostT)); % StimOn moment (0 for aud, 0.25 for som) - postT

        % convert to firing rate (Hz)
        baseFR(trial, cluster) = baseSpikes/abs(PreT);
        expFR(trial, cluster) = expSpikes/abs(PostT);
    end
end
end
