%% Quantify results
% Data paths to sorted data, behavioral data file, output folder
% To be done after PostCuration_BK.m
% Post processing data steps:
%   1. get firing rate during set time frame for each trial of each unit
%   2. quatify which cluster responds to which event
%       responsive unit = sig different firing after stim vs baseline
%   3. quantify response strength, make various plots
%   4. combine data from multiple animals

clearvars

%  set directories
session = 2;

%recordingFolder = 'D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\';
%acfeedbackPath = [recordingFolder 'Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B'];
%recPath = [recordingFolder 'Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\'];
%TTLPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\'];
%messagesPath = [recordingFolder 'Record Node 103\experiment1\recording1\events\MessageCenter\'];
%rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M8\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M8'; % output directory

Fs = 30000; % sampling freq

[cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, StimOn] = loadData(OutPath, session, BehaviorPath);

%% 1. calculate firing rates

% define pre and post simulus onset
PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/4)/1000; % baseline period
if strcmp(stimuli_parameters.Par.Rec, "SxA")
    PostT = 0.25; % captures initial noise period & half of vibrotac (in noise) period
elseif strcmp(stimuli_parameters.Par.Rec, "SOM") && strcmp(stimuli_parameters.Par.SomatosensoryWaveform, "Square")
    PostT = 0.05; % best way to capture onset stimulus
elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
    PostT = 0.2; % limited by previous recordings
else
    disp("No post stimulus time frame defined")
end

% calculate firing rate (Hz) in time window
[baselineRate, stimulusRate] = firingrate(aligned_spikes, PreT, PostT, StimOn);

% calculate stimulus induced change in firing rate
dstimulusRate = stimulusRate - baselineRate;
%% 2. quantify reactive units
% ! edit exp condition for SxA sessions
% ! edit AMn indexing if needed (needed in earlier data sets)

[responsive_units, resp_cids] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids);

%% 3. calculate mean stimulus induced change in firing

% select correct amplitude and frequency parameters
if strcmp(stimuli_parameters.Par.Rec, "AMn") % noise session
    uAmp = unique(stimuli_parameters.Stm.Intensity); 
    uFreq = unique(stimuli_parameters.Stm.Mf);
else %SxA and SOM sessions
    uAmp = unique(stimuli_parameters.Stm.Amplitude);
    uFreq = unique(stimuli_parameters.Stm.SomFreq);
end

nAmp = length(uAmp);
nFreq = length(uFreq);

% select multiple stimuli types in one session
if strcmp(stimuli_parameters.Par.Rec, "SxA") %dual modal session
    conditions = {"OO", "OA", "SO", "SA"};
else % unimodal sessions
    conditions = 1;
end
% initialize output variables
firing_mean = nan(length(conditions), nFreq, nAmp, length(resp_cids));
firing_se = nan(length(conditions), nFreq, nAmp, length(resp_cids));

% mean dfiring rate per stimulus combination for all responsive units
for condition = 1:length(conditions)
    for freq = 1:nFreq
        for amp = 1:nAmp
            index = stimuli_parameters.Stm.Amplitude == uAmp(amp) & stimuli_parameters.Stm.SomFreq == uFreq(freq) & condition;
            firing_mean(condition, freq, amp, :) = mean(dstimulusRate(index, resp_cids));
            firing_se(condition, freq, amp, :) = mean(dstimulusRate(index, resp_cids));
        end
    end
end

%% 3.5 Save if needed
filename = sprintf('M%.2i_S%.2i_%s_ResponseProperties', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
save(fullfile(OutPath, filename), "responsive_units", "baselineRate", "stimulusRate", "firing_mean", "firing_se")

%% 4. combine datasets across mice

% animals to include
animal = [4 6 7];

%define session per animal
session = 2;
%session_options = {"AMn", "SOM", "SxA", "FRA"};
%session = char(session_options{1});

% get data from session to analyse
for file = 1:length(animal)
    OutPath = ['D:\DATA\Processed\M' num2str(animal(file), '%d') '\'];
    sessionFile = ['M' num2str(animal(file), '%.2d') '_S' num2str(session, '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
    stim_files = dir(fullfile(OutPath, sessionFile));
    files = load([stim_files.folder '\' stim_files.name]);
    
    % concatinate data from all units
    data_all = [data_all; files.dfiring_mean'];
end

%% plot for multiple animals
uAmp = [0 15 30 45 60];

data = mean(data_all);
errors = std(data_all) / sqrt(size(data_all, 1));

figure;
hold on
plot(data, 'LineWidth', 1.25)
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(1:length(uAmp))
xticklabels(uAmp)

ylabel('df Firing rate (Hz)')
xlabel('Broadband noise intensity (dB SPL)')

%% quantify response strength per condition

conditions = {"OO", "OA", "SO", "SA"}; % defines order
firing_mean = nan(length(conditions), length(resp_cids));
firing_se= nan(length(conditions), length(resp_cids));

for i = 1:length(conditions)
    index = strcmp(stimuli_parameters.Stm.MMType, conditions{i});
    firing_mean(i, :) = mean(stimulusRate(index, resp_cids));
    firing_se(i,:) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
end

%firing_mean(condition, freq, amp, :)

%plot bar graph
data = mean(firing_mean, 2);
errors = mean(firing_se, 2);

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

%% pressure tuning curve

uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
nClusters = length(cids);
dfiring = nan(nAmp, nClusters);

figure;
hold on

for cluster = 1:nClusters

    for amp = 1:nAmp
        index = stimuli_parameters.Stm.Amplitude == uAmp(amp);
        dfiring(amp, cluster) = mean(dstimulusRate(index, cluster));
    end
    
    % normalize to control
    dfiring(2:5, cluster) = dfiring(2:5, cluster) - dfiring(1, cluster);

    % plot responsive units in distinct colour
    if any(ismember(resp_cids, cluster))
        plot(dfiring(:, cluster), 'r');
    else
        plot(dfiring(:, cluster), 'k');
    end

    ylabel('df Firing rate (Hz)')
    xlabel('stimulus strength (V)')
    xticks(1:nAmp)
    xticklabels(uAmp)
end
hold off

% pressure tuning curve mean of only responsive units
firing_mean = nan(length(conditions), length(resp_cids));
firing_se= nan(length(conditions), length(resp_cids));

for amp = 1:nAmp
    index = stimuli_parameters.Stm.Amplitude == uAmp(amp);
    firing_mean(amp, :) = mean(stimulusRate(index, resp_cids));
    firing_se(amp,:) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
end

% plot indivudual traces along average
% plot(dfiring_mean, 'r')
% plot(mean(dfiring_mean, 2), 'k', 'LineWidth', 1.25)
% xticks(1:nAmp)
% xticklabels(uAmp)
% 
% ylabel('df Firing rate (Hz)')
% xlabel('stimulus strength (V)')

figure;
hold on
data = mean(dfiring_mean, 2); %mean(firing_mean, 2);
errors = mean(dfiring_se, 2); %mean(firing_sd, 2);
plot(data, 'LineWidth', 1.25)
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(1:nAmp)
xticklabels(uAmp)

ylabel('df Firing rate (Hz)')
xlabel('stimulus strength (V)')

%% vibrotac tuning curves

uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
uFreq = unique(stimuli_parameters.Stm.SomFreq);
nFreq = length(uFreq);

nClusters = length(cids);

figure;
hold on

for cluster = 1:nClusters

    for freq = 1:nFreq
        for amp = 1:nAmp
            index = stimuli_parameters.Stm.Amplitude == uAmp(amp) & stimuli_parameters.Stm.SomFreq == uFreq(freq);
            dfiring(freq, amp, cluster) = mean(dstimulusRate(index, cluster));
        end
    end

    % normalize to control condition
    dfiring(2:nFreq, 2:nAmp, cluster) = dfiring(2:nFreq, 2:nAmp, cluster) - dfiring(1, 1, cluster);

    % plot responsive units in distinct colour
    if any(ismember(resp_cids, cluster))
        plot(dfiring(2:nFreq, 2:nAmp, cluster), 'r');
    else
        plot(dfiring(2:nFreq, 2:nAmp, cluster), 'k');
    end
end

ylabel('df Firing rate (Hz)')
xlabel('Vibrotactile frequency (Hz)')
xticklabels(uFreq(2:nFreq))
hold off

% vibrotac tuning curve mean of only responsive units
% works for SxA and SOM sessions
dfiring_mean = nan(nFreq, nAmp, resp_cids);
dfiring_se = nan(nFreq, nAmp, resp_cids);

for freq = 1:nFreq
    for amp = 1:nAmp
        index = stimuli_parameters.Stm.Amplitude == uAmp(amp) & stimuli_parameters.Stm.SomFreq == uFreq(freq);
        dfiring_mean(freq, amp, :) = mean(stimulusRate(index, resp_cids));
        dfiring_se(freq, amp, :) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
    end
end

% normalize to control condition
dfiring_mean(2:nFreq, 2:nAmp, :) = dfiring_mean(2:nFreq, 2:nAmp, :) - dfiring_mean(1, 1, :);

figure;
hold on
data = mean(dfiring_mean(2:nFreq, 2:nAmp, :), 3); %mean(firing_mean, 2);
errors = mean(dfiring_se(2:nFreq, 2:nAmp, :), 3); %mean(firing_sd, 2);
plot(data, 'LineWidth', 1.25)
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

xticklabels(uFreq(2:nFreq))
legend(num2str(uAmp(2:end)))
ylabel('df Firing rate (Hz)')
xlabel('Vibrotactile frequency (Hz)')


%% noise intensity level tuning curve
uAmp = unique(stimuli_parameters.Stm.Intensity);
nAmp = length(uAmp);
nClusters = length(cids);

% pressure tuning curve mean of only responsive units
for amp = 1:nAmp
    index = stimuli_parameters.Stm.Intensity == uAmp(amp);
    dfiring_mean(amp, :) = mean(stimulusRate(index, resp_cids));
    dfiring_se(amp,:) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
end

figure;
hold on
data = mean(dfiring_mean, 2); %mean(firing_mean, 2);
errors = mean(dfiring_se, 2); %mean(firing_sd, 2);
plot(data, 'LineWidth', 1.25)
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(1:nAmp)
xticklabels(uAmp)

ylabel('df Firing rate (Hz)')
xlabel('Broadband noise intensity (dB SPL)')

data = mean(data_all);
errors = std(data_all) / sqrt(size(data_all, 1));
figure;
hold on
plot(data, 'LineWidth', 1.25)
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(1:length(uAmp))
xticklabels(uAmp)

ylabel('df Firing rate (Hz)')
xlabel('Broadband noise intensity (dB SPL)')

%% compare firing rate unimodal to multimodal stimulus presentation
% makes scatter plot of firing rates per responsive unit between two conditions

conditions = {"OO", "OA", "SO", "SA"};
dfiring = nan(length(conditions), size(dstimulusRate, 2));

for i = 1:length(conditions)
    index = strcmp(stimuli_parameters.Stm.MMType, conditions{i});
    dfiring_mean(i, :) = mean(dstimulusRate(index, resp_cids));
    dfiring_se(i,:) = std(dstimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
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
index = strcmp(stimuli_parameters.Stm.MMType, "SO"); % to do: select best freq. use 3rd dimention matrix as condition index
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


%% best modality per unit
figure;
bar([length(audPref) length(somPref)])

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
%% ------------------------- Local functions ------------------------- %%
function [cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, StimOn] = loadData(OutPath, session, BehaviorPath)
% get data files
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
else
    disp("Extract Srise before continuing")
end

aligned_spikes = aligned_spikes.SpkT;

% load sessions details
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']).name;
sessions_TTLs = load([OutPath '\' TTLs_file]);

% add delay to "SO"
if strcmp(stimuli_parameters.Par.Rec, 'SxA')
    StimOn = stimuli_parameters.Stm.SomAudSOA./1000;
    StimOn(isnan(stimuli_parameters.Stm.SomAudSOA)) = 0;
else
    StimOn = zeros(size(stimuli_parameters.Stm, 1), 1);
end
end