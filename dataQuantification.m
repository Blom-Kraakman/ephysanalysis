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


%SxA, SOM
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
%session = [4, 6]; % M9


%  set directories
animal = 9;
session = 10;

BehaviorPath = ['D:\DATA\Behavioral Stimuli\M' num2str(animal)]; % stimuli parameters
OutPath = ['D:\DATA\Processed\M' num2str(animal)]; % output directory

[cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, StimOn] = loadData(OutPath, session, BehaviorPath);

% 1. calculate firing rates

% Initiate table
Cluster = cids';
MouseNum = repmat(animal,[length(Cluster),1]);
%Session = repmat(session,[length(Cluster),1]);
unitResponses = table(MouseNum, Cluster);

% define pre and post simulus onset
if strcmp(stimuli_parameters.Par.Rec, "SxA")
    %PostT = 0.25; % captures initial noise period & half of vibrotac (in noise) period
    PostT = 0.5; % whole vibrotac + dual mdoe period
    PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
    %StimOn = zeros(length(stimuli_parameters.Stm.SomAudSOA),1) + 0.250;
elseif strcmp(stimuli_parameters.Par.Rec, "SOM") && strcmp(stimuli_parameters.Par.SomatosensoryWaveform, "Square")
    PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
    PostT = 0.1; % best way to capture onset stimulus
elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
    PostT = 0.2; % limited by previous recordings
    %PreT = (str2double(stimuli_parameters.Par.AMISI)/3)/1000; % baseline period
    PreT = 0.2;
else
    error("No analysis window defined")
end

% calculate firing rate (Hz) in time window
[baselineRate, stimulusRate] = firingrate(aligned_spikes, PreT, PostT, StimOn);

% calculate stimulus induced change in firing rate
dstimulusRate = stimulusRate - baselineRate;

% 2. calculate mean stimulus induced change in firing

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
    conditions = {'OO', 'OA', 'SO', 'SA'};
else % unimodal sessions
    conditions = 1;
end

% initialize output variables
firing_mean = nan(nFreq, nAmp, length(conditions), length(cids));

% mean dfiring rate per stimulus combination for all units
for condition = 1:length(conditions)
    for freq = 1:nFreq
        for amp = 1:nAmp
            if strcmp(stimuli_parameters.Par.Rec, "SxA")
                index = stimuli_parameters.Stm.Amplitude == uAmp(amp) & stimuli_parameters.Stm.SomFreq == uFreq(freq) & strcmp(stimuli_parameters.Stm.MMType, conditions{condition});
            elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
                index = stimuli_parameters.Stm.Intensity == uAmp(amp) & stimuli_parameters.Stm.Mf == uFreq(freq);
            else
                index = stimuli_parameters.Stm.Amplitude == uAmp(amp) & stimuli_parameters.Stm.SomFreq == uFreq(freq);
            end

            firing_mean(freq, amp, condition, :) = mean(dstimulusRate(index, 1:length(cids)));

        end
    end
end
clear freq amp condition

% 3. quantify reactive units
% ! edit AMn indexing if needed (needed in earlier data sets)
for cond = 1:length(conditions)
    responsive = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, conditions(cond));
    
    % index responsive units
    idx = max(responsive == unitResponses.Cluster, [], 2);

    if isempty(idx)
        idx = false(length(cids),1);
    end

    % add responsive true/false per cluster per condition
    unitResponses = addvars(unitResponses, idx);

end

if strcmp(stimuli_parameters.Par.Rec, "SxA")
    % add OA sound steady state
    StimOn = zeros(length(stimuli_parameters.Stm.SomAudSOA),1) + 0.250;
    [baselineRate, stimulusRate] = firingrate(aligned_spikes, PreT, PostT, StimOn); % calculate firing rate (Hz) in time window
    dstimulusRate = stimulusRate - baselineRate; % calculate stimulus induced change in firing rate

    %calculate mean stimulus induced change in firing
    firing_mean = nan(nFreq, nAmp, length(conditions), length(cids));

    % mean dfiring rate per stimulus combination for all units
    for condition = 1:length(conditions)
        for freq = 1:nFreq
            for amp = 1:nAmp
                index = stimuli_parameters.Stm.Amplitude == uAmp(amp) & stimuli_parameters.Stm.SomFreq == uFreq(freq) & strcmp(stimuli_parameters.Stm.MMType, conditions{condition});
                firing_mean(freq, amp, condition, :) = mean(dstimulusRate(index, 1:length(cids)));
            end
        end
    end

    responsive = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, 'OA');

    % index responsive units
    idx = max(responsive == unitResponses.Cluster, [], 2);

    if isempty(idx)
        idx = false(length(cids),1);
    end

    % add responsive true/false per cluster per condition
    unitResponses = addvars(unitResponses, idx);

    % modify table comlumn names and add to struct

    unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster', 'OO', 'OA', 'SO', 'SA', 'OA+0.25'};

elseif strcmp(stimuli_parameters.Par.Rec, "SOM")
    unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster', 'SOM'};
end

StimResponseFiring.unitResponses = unitResponses;
StimResponseFiring.MouseNum = stimuli_parameters.Par.MouseNum;
StimResponseFiring.session = stimuli_parameters.Par.Set;
StimResponseFiring.type = stimuli_parameters.Par.Rec;
StimResponseFiring.conditions = conditions;
StimResponseFiring.cids = cids;
StimResponseFiring.ampitudes = uAmp;
StimResponseFiring.frequencies = uFreq;
StimResponseFiring.PreT = PreT;
StimResponseFiring.PostT = PostT;
StimResponseFiring.StimOn = StimOn;
StimResponseFiring.baselineRate = baselineRate;
StimResponseFiring.stimulusRate = stimulusRate;
StimResponseFiring.dfiringRate = dstimulusRate;
StimResponseFiring.firing_mean = firing_mean;

% 3.5 Save if needed
filename = sprintf('M%.2i_S%.2i_%s_ResponseProperties', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
save(fullfile(OutPath, filename), "StimResponseFiring")

%% 4. add SOM responses to table
animal = 9; % M11
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
session = [4, 6]; % SxA, SOM M9

OutPath = ['D:\DATA\Processed\M' num2str(animal, '%d') '\ResponseProperties'];

% load SxA data
sessionFile = ['M' num2str(animal, '%.2d') '_S' num2str(session(1), '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
stim_files = dir(fullfile(OutPath, sessionFile));
dataS = load([stim_files.folder '\' stim_files.name]);
unitResponses_SxA = dataS.StimResponseFiring.unitResponses;

% load pressure data
sessionFile = ['M' num2str(animal, '%.2d') '_S' num2str(session(2), '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
stim_files = dir(fullfile(OutPath, sessionFile));
dataS = load([stim_files.folder '\' stim_files.name]);
unitResponses_SOM = dataS.StimResponseFiring.unitResponses;

% combine tables
unitResponses = joindata(unitResponses_SxA, unitResponses_SOM);

%include condition when saving
filename = sprintf('M%.2i_S%.2i&S%.2i_UnitResponses', animal, session(1), session(2));
save(fullfile(OutPath, filename), "unitResponses")

%% 5. combine datasets across mice

% animals to include
%animal = [4 6 7];
animal = [9 10 11];

%define session per animal
%session = [2 2 2]; % noise M4,6,7
session = [10 6 6]; % SxA earplug
%session = [6 5 4]; % step M9-11
%session = [5 4 2 3]; %SxA M8-11
data_all = [];
unitResponses_all = table;

% get data from session to analyse
for file = 1:length(animal)

    OutPath = ['D:\DATA\Processed\M' num2str(animal(file), '%d')];

    %OutPath = ['D:\DATA\Processed\M' num2str(animal(file), '%d') '\ResponseProperties\'];
    sessionFile = ['M' num2str(animal(file), '%.2d') '_S' num2str(session(file), '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
    stim_files = dir(fullfile(OutPath, sessionFile));
    dataS = load([stim_files.folder '\' stim_files.name]);

    %disp(['post stim time: ' num2str(dataS.StimResponseFiring.PostT) 'sec'])
    %disp(['pre stim time: ' num2str(dataS.StimResponseFiring.PreT) 'sec'])

    % concatinate data from all units
    % select common amp and freq data points
    if strcmp(dataS.StimResponseFiring.type, 'SxA')

        if animal(file) == 9 && session(file) == 4
            data = dataS.StimResponseFiring.firing_mean(:, [1,3,5], :, :);
        elseif animal(file) == 8 || animal(file) == 10 || animal(file) == 11
            data = dataS.StimResponseFiring.firing_mean(1:8, :, :, :);
        else
            data = dataS.StimResponseFiring.firing_mean; % no changes needed
        end
    else
        data = dataS.StimResponseFiring.firing_mean; % no changes needed
    end

    data_all = cat(4, data_all, data);

    % combine unitResponses tables
    if strcmp(dataS.StimResponseFiring.type, 'SxA')
        sessionFile = ['M' num2str(animal(file), '%.2d') '_S' num2str(session(file), '%.2d') '*UnitResponses.mat']; % select based on stimulus type '_*.mat'
        stim_files = dir(fullfile(OutPath, sessionFile));
        unitResponses = load([stim_files.folder '\' stim_files.name]);
        unitResponses_all = vertcat(unitResponses_all, unitResponses.unitResponses);
    else
        unitResponses_all = vertcat(unitResponses_all, dataS.StimResponseFiring.unitResponses);
    end

end

%% save if needed
filename = sprintf('M09-10-11_SxA_earplug_UnitResponses');
save(fullfile(OutPath, filename), "unitResponses_all", "data_all")

%% PLOTTING START


%% combine unit responses
unitResponses = unitResponses_all;
sound = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & ~unitResponses.SO & ~unitResponses.SOM);
vibrotac = table2array(unitResponses.SO & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SOM);
pressure = table2array(unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO);
sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
sound_pressure = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & ~unitResponses.SO);
vibrotac_pressure = table2array((unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & unitResponses.SO));
all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);
non = table2array(~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SA & ~unitResponses.SO);

clear data dataS unitResponses_all
%% venn diagram response groups
sets = ["sound" "vibrotactile" "pressure"];
labels = [sum(sound), sum(vibrotac), sum(pressure), sum(sound_vibrotac), sum(sound_pressure), sum(vibrotac_pressure), sum(all)]; % order: A, B, C, A&B, A&C, B&C and A&B&C

venn(3,'sets',sets,'labels',labels, 'colors', [0.247 0.635 0.831; 0.266 0.666 0.6;0.808 0.529 0.666], 'alpha',0.75);

%% quantify response strength per condition
% data format:
%   firing_mean(freq, amp, condition, :)
   conditions = ["OO", "OA", "SO", "SA"];

% select all units responding to sound and vibration
unitResponses = unitResponses_all;
sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
%sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
%all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);

% select and format data
sound_vibro_idx = max(all, sound_vibrotac);
%data = (mean(data, 1, "omitnan"));
%data = squeeze(mean(mean(data, 1, "omitnan"), 2, "omitnan"));

% make data matrix to plot
data_OO = squeeze(data_all(1,1,1,sound_vibro_idx));
data_OA = squeeze(data_all(1,1,2,sound_vibro_idx));
%data_SO = squeeze(data_all(3,3,3,sound_vibro_idx)); % 20Hz, 0.3V
%data_SA = squeeze(data_all(3,3,4,sound_vibro_idx));

[maxValue, maxValueIdx] = max(max(data_all(:,:,3,sound_vibro_idx),[],2));
maxValueIdx = squeeze(maxValueIdx);
sound_vibro_units = find(sound_vibro_idx);
data_SO = zeros(length(maxValueIdx),1);
data_SA = zeros(length(maxValueIdx),1);
for i = 1:length(maxValueIdx)
    data_SO(i,1) = squeeze(data_all(maxValueIdx(i),3,3,sound_vibro_units(i))); % max Hz, 0.3V
    data_SA(i,1) = squeeze(data_all(maxValueIdx(i),3,4,sound_vibro_units(i)));
end

data = [data_OO'; data_OA'; data_SO'; data_SA']; % 20Hz, 0.3V
data(5,:) = data(2,:) + data(3,:);
%data = log(data);
%errors = std(data(:, :), 0, 2); % / sqrt(size(data, 2));

%plot bar graph
figure;
set(gcf,'position',[500,150,900,700])
hold on
fig = bar(mean(data,2)); % avg FR of all responsive units
x = repmat((1:5)',1,length(maxValueIdx));
for i = 1:5
    swarmchart(x(i,:), data(i,:), 20, 'k', 'filled','XJitterWidth',0.4);
end

%plot(categorical(conditions),data, '.k', 'MarkerSize', 15)
%plot(1:4,data,'o', 'Color', [0.5 0.5 0.5])
%errorbar(1:size(data), mean(data,2), errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
%yline(data(4))
%scatter(1:length(data), firing_mean, 10, 'k', 'filled');

% bar colors
set(fig, 'FaceColor', 'Flat')
fig.CData = [0 0 0; 0.808 0.529 0.666; 0.247 0.635 0.831; 0.580 0.455 0.651; 0.5 0.5 0.5];

% axis
set(gca,'fontsize',18)
%xlabel('conditions')
xticks(1:5)
xticklabels(["control", "sound", "vibrotactile", "multimodal", "sound + vibrotactile"])
ylabel('\Delta Firing rate (spikes/s)')
%title(PostT)

hold off

%% calculate preference and modulation index 

% select all units responding to sound and vibration
unitResponses = unitResponses_all;
sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
%sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
%all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);
%non = table2array(~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SA & ~unitResponses.SO);
non = table2array(~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SA & ~unitResponses.SO);

% select and format data
sound_vibro_idx = max(all, sound_vibrotac);
%data = data_all(:,:,:,sound_vibro_idx);
%aud_data = abs(squeeze(data(1,1,2,:)));
%som_data = abs(squeeze(mean(data(2:7,3,3,:), 1)));
data = data_all(:,:,:,~non);
data = squeeze(mean(mean(data, 1, "omitnan"), 2, "omitnan"));
%data(5,:) = abs(data(2,:)) + abs(data(3,:));

pref_index = (data(2,:) - data(3,:)) ./ (data(2,:) + data(3,:));

pref_index(isnan(pref_index)) = [];
pref_index(pref_index > 80) = [];
mean(pref_index)
std(pref_index)


[N, edges] = histcounts(pref_index, -6:0.5:2);
%[N, edges] = histcounts(modulation_index);
figure;
for i = 1:(length(edges)-1)
    binCenters(:,i) = mean(edges(:,i:i+1));
    i = i+1;
end

bar(binCenters, N, 1, 'FaceColor', [0.5 0.5 0.5])
xline(0, 'k')
xlabel('Modality preference')
ylabel('# units')
set(gca,'fontsize',18)
clear binCenters

modulation_index = data(4,:) - (data(2,:) + data(3,:));
modulation_index(isnan(modulation_index)) = 0;

mean(modulation_index)
std(modulation_index)

% plot(sort(modulation_index))
%yline(0)
%xlabel('unit');
%ylabel('modualtion index')

[N, edges] = histcounts(modulation_index, -3:0.5:6);
%[N, edges] = histcounts(modulation_index);
figure;
for j = 1:(length(edges)-1)
    binCenters(:,j) = mean(edges(:,j:j+1));
    j = j+1;
end

bar(binCenters, N, 1, 'FaceColor', [0.5 0.5 0.5])
xline(0, 'k')
xlabel('Modulation strength')
ylabel('# units')
set(gca,'fontsize',18)

clear binCenters

%% compare firing rate unimodal to multimodal stimulus presentation
% makes scatter plot of firing rates per responsive unit between two conditions
% data format:
%   firing_mean(freq, amp, condition, resp_cids)
%   conditions = {"OO", "OA", "SO", "SA"};

% select all units responding to sound and vibration
unitResponses = unitResponses_all;
sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
%sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
%all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);

% select and format data
sound_vibro_idx = max(all, sound_vibrotac);
data = data_all(:,:,:,sound_vibro_idx);

%data = data_all;
%data = squeeze(mean(mean(data, 1, "omitnan"), 2, "omitnan"));
aud_data = abs(squeeze(data(1,1,2,:)));
som_data = abs(squeeze(mean(data(2:7,3,3,:), 1)));

% for unit = 1:size(data,4)
% 
%     % stat test
%     [p,h,stats] = signrank(aud_data, som_data, 'alpha', 0.01); % different from control trials?
% 
%     if h
%         if abs(aud_trials(unit)) > abs(som_trials(unit))
%             audPref = [audPref; unitResponses.responsive(unit)];
%         elseif abs(aud_trials(unit)) < abs(som_trials(unit))
%             somPref = [somPref; unitResponses.responsive(unit)];
%         elseif abs(aud_trials(unit)) == abs(som_trials(unit))
%             noPref = [noPref; unitResponses.responsive(unit)];
%         end
%     else
%         noPref = [noPref; unitResponses.responsive(unit)];
%     end
% end

% best modality per unit
figure;
hold on
scatter(aud_data, som_data, 'k', "filled")
% add 45 degree angle line
plot([-1, 15],[ -1, 15], 'k')

% format
xlabel('\Delta Firing rate sound trials (spikes/s)')
ylabel('\Delta Firing rate vibrotactile trials (spikes/s)')
axis([-1 15 -1 15])

set(gca,'fontsize',18)

%% plot noise for multiple animals
uAmp = [0 15 30 45 60];

data_all = squeeze(data_all);
data = mean(data_all,2);
%errors = std(data_all,[],2) / sqrt(size(data_all, 2));

figure;
hold on
plot(uAmp, data_all, 'Color', "#CE87AA")
plot(uAmp, data, 'Color', "#9E0E56", 'LineWidth', 3)

%errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(uAmp)

ylabel('\Delta Firing rate (Hz)')
xlabel('Broadband noise intensity (dB SPL)')
set(gca,'fontsize',18)

%% pressure tuning curve
% work with 4D data format? 

% variables
uAmp = unique(StimResponseFiring.ampitudes)';
nAmp = length(uAmp);

umN = round((uAmp * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

% format data
data = squeeze(data_all)';
index = unitResponses_all.SOM;

FRpressure_mean = mean(data(index,:));
FRpressure_mean_c = mean(data(~index,:));

FRpressure_med = median(data(index,:));
FRpressure_med_c = median(data(~index,:));

%FRpressure_std = std(data(index,:));
%FRpressure_std_c = std(data(~index,:));

FRpressure_se = std(data(index,:)) / sqrt(size(data(index,:), 1));
FRpressure_se_c = std(data(~index,:)) / sqrt(size(data(~index,:), 1));

FRpressure_iqr = iqr(data(index,:));
FRpressure_iqr_c = iqr(data(~index,:));

%plot figure
figure;
hold on
plot(umN, data(index,:), 'Color',  "#44AA99", 'LineWidth', 0.1)
plot(umN, FRpressure_mean, 'Color',  "#267165", 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med, FRpressure_iqr, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

plot(umN, data(~index,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
plot(umN, FRpressure_mean_c, 'k', 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med_c, FRpressure_iqr_c, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% format figure
xlim([0 umN(nmN)])
set(gca,'fontsize',16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Stimulus strength (mN)')
%xticks(1:nAmp)
%xticklabels(uAmp)

%% vibrotac tuning curves
% of only responsive units
% work with 4D data format?

%uAmp = unique(stimuli_parameters.Stm.Amplitude);
uAmp = [0;0.1;0.3];
nAmp = length(uAmp);

umN = round((uAmp * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

%uFreq = unique(stimuli_parameters.Stm.SomFreq);
uFreq = [10;20;50;100;200;300;400];
nFreq = length(uFreq);

% select correct condition from SxA session
condition = 3; % ["OO", "OA", "SO", "SA"]

%responsive units
units_idx = unitResponses_all.SO;
data = data_all(2:8, 2:3, condition, units_idx);
FRvibrotac_mean = mean(data, 4, "omitnan");
FRvibrotac_se = std(data, [],4) / sqrt(size(data, 4));
FRvibrotac_std = std(data, [],4);

%non responsive units
no_resp_idx = ~unitResponses_all.SO;
data_c = data_all(2:8, 2:3, condition, no_resp_idx);
FRvibrotac_mean_c = mean(data_c, 4, "omitnan");
FRvibrotac_se_c = std(data_c, [],4) / sqrt(size(data_c, 4));
FRvibrotac_std_c = std(data_c, [],4);

% plot
figure;
hold on

%errorbar(uFreq, FRvibrotac_mean, FRvibrotac_se, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
for cluster = 1:size(data,4)
    %plot(uFreq , data(:, 1, cluster), 'Color', [0.5, 0.5, 0.5])
    plot(uFreq , data(:, 2, cluster), 'Color', "#8BC9E8")

end

for cluster_c = 1:size(data_c,4)
    plot(uFreq , data_c(:, 2, cluster_c), 'Color', [0.5, 0.5, 0.5])
end

plot(uFreq, FRvibrotac_mean(:, 2), 'LineWidth', 3, 'Color', "#3FA2D4")
%errorbar(uFreq, FRvibrotac_mean, FRvibrotac_std, 'b', 'linestyle', 'none', 'LineWidth', 0.5);

plot(uFreq, FRvibrotac_mean_c(:, 2), 'k', 'LineWidth', 3)
%errorbar(uFreq, FRvibrotac_mean_c, FRvibrotac_std_c, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% format figure
%plot(FRvibrotac_mean(2:nFreq, 2:nAmp));

ax = gca;
xlim = ([10 400]);
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Vibrotactile frequency (Hz)')
%legend(num2str(umN(2:end)))
set(ax,'fontsize',18)
xticks([0 10 20 50 100 200 300 400])
ax.XScale = 'log';


%% vibrotac tuning curve single units
% of only responsive units
% work with 4D data format?

%uAmp = unique(stimuli_parameters.Stm.Amplitude);
uAmp = [0;0.1;0.3];
nAmp = length(uAmp);
umN = round((uAmp * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

%uFreq = unique(stimuli_parameters.Stm.SomFreq);
uFreq = [10;20;50;100;200;300;400];
nFreq = length(uFreq);

% select correct condition from SxA session
condition = 3; % ["OO", "OA", "SO", "SA"]

%responsive units
%units_idx = unitResponses_all.SO;
unit_idx = 41;
data = data_all(2:8, 2:3, condition, unit_idx);

% plot
figure;
hold on

plot(uFreq , data(:, 1, cluster), 'LineWidth', 3, 'Color', "#8BC9E8")
plot(uFreq , data(:, 2, cluster), 'LineWidth', 3, 'Color', "#3FA2D4")

ax = gca;
%xlim([10 400]);
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Vibrotactile frequency (Hz)')
legend('16mN', '47mN')
set(ax,'fontsize',18)
xticks([0 10 20 50 100 200 400])
ax.XScale = 'log';

%% pressure tuning curve single unit 

% variables
uAmp = unique(StimResponseFiring.ampitudes)';
nAmp = length(uAmp);
umN = round((uAmp * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

% format data
data = squeeze(data_all)';
index = unitResponses_all.SOM;

FRpressure_mean = mean(data(index,:));
FRpressure_mean_c = mean(data(~index,:));

FRpressure_med = median(data(index,:));
FRpressure_med_c = median(data(~index,:));

%FRpressure_std = std(data(index,:));
%FRpressure_std_c = std(data(~index,:));

FRpressure_se = std(data(index,:)) / sqrt(size(data(index,:), 1));
FRpressure_se_c = std(data(~index,:)) / sqrt(size(data(~index,:), 1));

FRpressure_iqr = iqr(data(index,:));
FRpressure_iqr_c = iqr(data(~index,:));

%plot figure
figure;
hold on
plot(umN, data(index,:), 'Color',  "#44AA99", 'LineWidth', 0.1)
plot(umN, FRpressure_mean, 'Color',  "#267165", 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med, FRpressure_iqr, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

plot(umN, data(~index,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
plot(umN, FRpressure_mean_c, 'k', 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med_c, FRpressure_iqr_c, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% format figure
xlim([0 umN(nmN)])
set(gca,'fontsize',16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Stimulus strength (mN)')
%xticks(1:nAmp)
%xticklabels(uAmp)

% %% noise intensity level tuning curve
% uAmp = unique(Stm.Intensity);
% nAmp = length(uAmp);
% cids = cpos.id;
% nClusters = length(cids);
% 
% % pressure tuning curve mean of only responsive units
% for amp = 1:nAmp
%     index = Stm.Intensity == uAmp(amp);
%     dfiring_mean(amp, :) = mean(stimulusRate(index, resp_cids));
%     dfiring_se(amp,:) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
% end
% 
% figure;
% hold on
% data = mean(dfiring_mean, 2); %mean(firing_mean, 2);
% errors = mean(dfiring_se, 2); %mean(firing_sd, 2);
% plot(data, 'LineWidth', 1.25)
% errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
% xticks(1:nAmp)
% xticklabels(uAmp)
% 
% ylabel('df Firing rate (spikes/s)')
% xlabel('Broadband noise intensity (dB SPL)')
% 
% data = mean(data_all);
% errors = std(data_all) / sqrt(size(data_all, 1));
% figure;
% hold on
% plot(data, 'LineWidth', 1.25)
% errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
% xticks(1:length(uAmp))
% xticklabels(uAmp)
% 
% ylabel('df Firing rate (spikes/s)')
% xlabel('Broadband noise intensity (dB SPL)')
% set(gca,'fontsize',14)
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

Fs = 30000;
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

% load cluster unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']); %.name;
if ~isempty(cpos_file)
    cpos = load([OutPath '\' cpos_file.name]);
    cids = cpos.cpos.id';
else
    cids = [];
    disp('no cluster details found')
end

% load corresponsing  behaviour files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

% load aligned spike times
aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
if ~isempty(aligned_spikes_files)
    aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
    
    if isfield(aligned_spikes,"Srise")
        disp("Srise/Sfall loaded from data file.")
        Srise = aligned_spikes.Srise;
        Sfall = aligned_spikes.Sfall;
    else
        Srise = [];
        Sfall = [];
        disp("Extract Srise before continuing")
    end

    aligned_spikes = aligned_spikes.SpkT;

else
    Srise = [];
    Sfall = [];
    aligned_spikes = [];
    disp('no aligned spikes file found')
end

% load sessions details
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']); %.name;
if ~isempty(TTLs_file)
    sessions_TTLs = load([OutPath '\' TTLs_file.name]);
else
    sessions_TTLs = [];
    disp('no session TTL file found')
end

% add delay to "SO"
if strcmp(stimuli_parameters.Par.Rec, 'SxA')
    StimOn = stimuli_parameters.Stm.SomAudSOA./1000;
    StimOn(isnan(stimuli_parameters.Stm.SomAudSOA)) = 0;
else
    StimOn = zeros(size(stimuli_parameters.Stm, 1), 1);
end
end