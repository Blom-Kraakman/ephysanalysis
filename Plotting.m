%% ---------------------- Load data  ---------------------- %%
% aligned spike data

close all
clearvars

OutPath = 'D:\DATA\Processed\M20\ICX'; % output directory
KSPath = 'D:\DATA\EphysRecordingsSorted\M20\ICX\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M20\'; % stimuli parameters

% load stimuli order files
session = 2;
[cids, stimuli_parameters, aligned_spikes, ~, ~, sessions_TTLs, onsetDelay, ~] = loadData(OutPath, session, BehaviorPath);

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

%% ---------------------- Channel map plotting  ---------------------- %%

% matches unit position to channel map
[~, channelpos] = plot_in_channel_map(KSPath, cpos.clusterinfo);

%% ---------------------- PostCuration Plotting ---------------------- %%
% plot aligned spikes from PostCuration_BK output files


%% plotting single sessions - in use
% raster plot of single units
% saves figures to OutPath

% select which session to plot
%session = 18;
for session = 13%relevant_sessions(1):relevant_sessions(2)
    % load corresponsing files
    sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
    stim_files = dir(fullfile(BehaviorPath, sessionFile));
    aligned_spikes_files = dir(fullfile(OutPath, sessionFile));

    if isempty(stim_files) || isempty(aligned_spikes_files)
        continue
    end

    stimuli_parameters = load([stim_files.folder '\' stim_files.name]);
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

    close all
end

%% plot single unit - in use
% raster plot + PSTH

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


%% plotting SOM sessions - keep
%saveplots = 0; %0 don't save, 1 save plots in OutPath
% currently plots raster + psth of each vibration freq seperatly

close all
%OutPath = 'D:\DATA\Processed\M10'; % output directory

% load unit info
%cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
%cpos = load([OutPath '\' cpos_file]);
%cids = cpos.cpos.id';

%sessions = 2; % [exp ctrl]
%load correct files
%session = sessions(1); % select experimental session

%sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
%stim_files = dir(fullfile(BehaviorPath, sessionFile));
%stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

%aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
%aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);

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

%% PSTH heatmap of all units - vibrotactile

% per stimulus combination, for all units
preT  = -0.2;
postT = 1.2;
binsize = 0.01; % 10ms
conditions = {'OO', 'SA', 'SO', 'OA'};

% set variables
amp = 0.3; % [0 0.1 0.3]
freq = 50;
%condition = 2; % 1: OO, 2: SA, 3: SO, 4: OA

% calculate PSTH per condition
for condition = 1:length(conditions)

    if condition == 1 || condition == 4
        dataidx = (stimuli_parameters.Stm.Var25 == condition);
    else
        dataidx = (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.SomFreq == freq) & (stimuli_parameters.Stm.Var25 == condition);
    end

    %for each cluster: spike count per bin
    clear spikecounts bin_edges
    for cluster = 1:length(cids)
        [N, edges] = histcounts(vertcat(aligned_spikes{dataidx, cluster}), preT:binsize:postT);
        spikecounts(cluster, :) = N;
        bin_edges(cluster, :) = edges;
    end
    bin_center = edges(1:end-1)+0.5*binsize;

    % baseline subtraction on spikecounts
    baseline = spikecounts(:,bin_center<0); %200ms baseline
    spikecounts = spikecounts - mean(baseline,2);

    % normalize spikecounts into PSTH
    % PSTH = sum spike counts / (number events * bin size)
    PSTH(:,:,condition) = spikecounts / (size(aligned_spikes(dataidx, cluster), 1) * binsize);

    % sort PSTH based on strongest response
    stimonset = 0;%find(edges == 0);
    stimoffset = 0.5;%find(edges == 0.5);
    [~, I] = sort(mean(PSTH(:,bin_center > stimonset & bin_center < stimoffset, condition),2), 1); % avg stimperiod firing
    PSTH(:,:,condition) = PSTH(I,:, condition);
end

% plot sound subtracted PSTH
figure;
subplot(3,1,1); % mulitmodal
plot(bin_center, mean(PSTH(:,:,2),1))
title('Condition: SA')
sgtitle([num2str(freq) 'Hz, ' num2str(amp) 'V'])
subplot(3,1,2); % sound only
plot(bin_center, mean(PSTH(:,:,4),1))
title('Condition: OA')
subplot(3,1,3); % avg PSTH of multimodal - sound only trace
plot(bin_center,(mean(PSTH(:,:,2),1)- mean(PSTH(:,:,4),1)))
title('Residual SA - OA')

condition = 3;
% plot individual + mean PSTH traces
figure;
subplot(2,1,1);
plot(bin_center,PSTH(:,:,condition))
title('individual traces')
subplot(2,1,2);
plot(bin_center,mean(PSTH(:,:,condition),1))
title('mean PSTH')
sgtitle(['Condition: ' conditions{condition} ' (' num2str(freq) 'Hz, ' num2str(amp) 'V)'])

%plot heatmap
figure;
imagesc(edges, (1:length(cids)), PSTH(:,:,condition))
colormap('hot'); % Choose a color map
cb = colorbar(gca, 'eastoutside');
cb.Label.String = 'spike/sec';
%clim([-5 30])
xlabel('Time (s)');
ylabel('clusters')
yticks(1:1:length(cids))
yticklabels(cids(I))
sgtitle(['Condition: ' conditions{condition} ' (' num2str(freq) 'Hz, ' num2str(amp) 'V)'])

%% PSTH heatmap of all units - Pressure

% per stimulus combination, for all units
preT  = -0.2;
postT = 0.5;
binsize = 0.01; % 10ms
conditions = {'OO', 'SA', 'SO', 'OA'};

% set variables
amp = 0.3160; % [0 0.013 0.032 0.063 0.095 0.127 0.19 0.253 0.316]
audint = 45;
%condition = 2; % 1: OO, 2: SA, 3: SO, 4: OA

% calculate PSTH per condition
for condition = 1:length(conditions)

    if condition == 1 || condition == 4
        dataidx = (stimuli_parameters.Stm.Var25 == condition);
    elseif condition == 3
        dataidx = (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.AudDur == 0) & (stimuli_parameters.Stm.Var25 == condition);
    else
        dataidx = (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.AudIntensity == audint) & (stimuli_parameters.Stm.Var25 == condition);
    end

    if ~sum(dataidx)
        error('no data selected')
    end

    %for each cluster: spike count per bin
    clear spikecounts bin_edges
    for cluster = 1:length(cids)
        [N, edges] = histcounts(vertcat(aligned_spikes{dataidx, cluster}), preT:binsize:postT);
        spikecounts(cluster, :) = N;
        bin_edges(cluster, :) = edges;
    end
    bin_center = edges(1:end-1)+0.5*binsize;

    % baseline subtraction on spikecounts
    baseline = spikecounts(:,bin_center<0); %200ms baseline
    spikecounts = spikecounts - mean(baseline,2);

    % normalize spikecounts into PSTH
    % PSTH = sum spike counts / (number events * bin size)
    PSTH(:,:,condition) = spikecounts / (size(aligned_spikes(dataidx, cluster), 1) * binsize);

    % sort PSTH based on strongest response
    stimonset = 0;%find(edges == 0);
    stimoffset = 0.1;%find(edges == 0.5);
    [~, I] = sort(mean(PSTH(:,bin_center > stimonset & bin_center < stimoffset, condition),2), 1); % avg stimperiod firing
    PSTH(:,:,condition) = PSTH(I,:, condition);
end

% plot sound subtracted PSTH
figure;
subplot(3,1,1); % mulitmodal
plot(bin_center, mean(PSTH(:,:,2),1))
title('Condition: SA')
sgtitle([num2str(amp) 'V'])
subplot(3,1,2); % sound only
plot(bin_center, mean(PSTH(:,:,4),1))
title('Condition: OA')
subplot(3,1,3); % avg PSTH of multimodal - sound only trace
plot(bin_center,(mean(PSTH(:,:,2),1)- mean(PSTH(:,:,4),1)))
title('Residual SA - OA')

condition = 3;
% plot individual + mean PSTH traces
figure;
subplot(2,1,1);
plot(bin_center,PSTH(:,:,condition))
title('individual traces')
subplot(2,1,2);
plot(bin_center,mean(PSTH(:,:,condition),1))
title('mean PSTH')
sgtitle(['Condition: ' conditions{condition} ' (' num2str(amp) 'V)'])

%plot heatmap
figure;
imagesc(edges, (1:length(cids)), PSTH(:,:,condition))
colormap('hot'); % Choose a color map
cb = colorbar(gca, 'eastoutside');
cb.Label.String = 'spike/sec';
%clim([-5 30])
xlabel('Time (s)');
ylabel('clusters')
yticks(1:1:length(cids))
yticklabels(cids(I))
sgtitle(['Condition: ' conditions{condition} ' (' num2str(amp) 'V)'])

%% ------------------- dataQuantification Plotting ------------------- %%
% plot aligned spikes from PostCuration_BK output files

OutPath = 'D:\DATA\Processed\M20\ICX';
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M20';
session = 2;
[~, ~, ~, ~, ~, ~, ~, StimResponseFiring] = loadData(OutPath, session, BehaviorPath);

%% -------------------------- PLOTTING START -------------------------- %%
%% combine unit responses
unitResponses = unitResponses_all;
sound = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & ~unitResponses.SO & ~unitResponses.SOM);
vibrotac = table2array(unitResponses.SO & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SOM);
pressure = table2array(unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO);
sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
sound_pressure = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & ~unitResponses.SO);
vibrotac_pressure = table2array((unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & unitResponses.SO));
all_stim = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);
non = table2array(~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SA & ~unitResponses.SO);

clear data dataS unitResponses_all

%% data selection new
stimOrder = readtable('D:\DATA\Processed\M12-16_stimOrder.csv');

%% venn diagram: response ratio groups
sets = ["sound" "vibrotactile" "pressure"];
total = size(data_all,4);
labels = [round((sum(sound)/total)*100), round((sum(vibrotac)/total)*100), ...
    round((sum(pressure)/total)*100), round((sum(sound_vibrotac)/total)*100),...
    round((sum(sound_pressure)/total)*100), round((sum(vibrotac_pressure)/total)*100), round((sum(all_stim)/total)*100)]; % order: A, B, C, A&B, A&C, B&C and A&B&C

venn(3,'sets',sets,'labels',labels, 'colors', [0.247 0.635 0.831; 0.266 0.666 0.6;0.808 0.529 0.666], 'alpha',0.75);

%% barplot: average response strength per condition
% data format:
%   firing_mean(freq, amp, condition, :)

conditions = ["OO", "OA", "SO", "SA"];
unitResponses = StimResponseFiring_all.unitResponses;
data_all = StimResponseFiring_all.firing_mean;

% select all units responding to sound and vibration
sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
%sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
all_stim = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
%all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);

% select and format data
sound_vibro_idx = max(all_stim, sound_vibrotac);
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
ylim([-40 80])
%title(PostT)

hold off

%% heatmap: calculate multisensory integration indexes 

% basic formulas:
% preference_index = (dFRsom - dFRaud) / (dFRsom + dFRaud)
% additivity_index = (dFRmulti - dFRsom - dFRaud) / (dFRsom + dFRaud)
% enhancement_magnitude = (dFRmulti - dFRmax) / dFRmax
% modulation_index = dFRmulti - (dFRsound + dFRsom)

data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids
signvalue = sign(data); % get sign of dFR
signvalue(signvalue(:,:,:) == 0) = 1; % make sign 1 if dFR=0

for c = 1:size(data,3)
    for i = 2:size(data,1)
        for j = 2:size(data,2)

            dFRsound = data(i,1,c);
            dFRvibrotac = data(1,j,c);
            dFRmulti = data(i,j,c);
            [dFRmax, max_index] = max([dFRsound, dFRvibrotac],[],"ComparisonMethod","abs");

            dFRmono_sign = [signvalue(i,1,c), signvalue(1,j,c)];
            dFRmulti_sign = signvalue(i,j,c);
            %dFRmax = max(abs(data(1,j,c)), abs(data(i,1,c)));
            %dFRmax = max(data(1,j,c),data(i,1,c),"ComparisonMethod","abs");

            % calculate multimodal enhancement index for each stim combination
            %enhancement_magnitude(i,j,c) = (dFRmulti - dFRmax) / dFRmax; % multisensory enhancement index
            enhancement_magnitude(i,j,c) = (((dFRmulti_sign * dFRmono_sign(max_index)) * abs(dFRmulti)) - abs(dFRmax)) / abs(dFRmax); % multisensory enhancement index = (dFRmulti - dFRmax) / dFRmax
            modulation_index(i,j,c) = dFRmulti - (dFRsound + dFRvibrotac); % multisensory modulation index, sign indicated direction multi
            additivity_index(i,j,c) = (dFRmulti - dFRvibrotac - dFRsound) / abs(dFRvibrotac + dFRsound);
            preference_index(i,j,c) = (dFRvibrotac - dFRsound) / (abs(dFRvibrotac) + abs(dFRsound)); % - pref sound; + pref vibrotac

        end
    end
end

clear i j c
% Plotting multisensory integration measures

dataToPlot = enhancement_magnitude;
%dataToPlot = modulation_index;
NaNindex_enhancement = ~isnan(enhancement_magnitude) & ~isinf(enhancement_magnitude);
NaNindex_modulation = ~isnan(modulation_index) & ~isinf(modulation_index);
NaNindex_additivity = ~isnan(additivity_index) & ~isinf(additivity_index);
NaNindex_preference = ~isnan(preference_index) & ~isinf(preference_index);

% plot enhancement magnitude index in heatmap
% variables
uparamA = unique(StimResponseFiring.amplitudes)';
%nAmp = length(uAmp);

umN = round((uparamA * 0.158)*1000); %0.158N/V callibation

for cluster = 12 %1:length(cids)
    figure('Position',[100,100,1000,800]);
    sgtitle(['unit ' num2str(cids(cluster))])

    % enhancement index
    subplot(2,2,1);
    imagesc(enhancement_magnitude(:,:, cluster),'AlphaData', NaNindex_enhancement(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    %clim([minV, maxV]);
    cb.Label.String = 'multisensory enhancement index';

    set(gca, 'YDir', 'normal', 'XTick',1:size(enhancement_magnitude,2),'XTickLabel',umN,'YTick',1:size(enhancement_magnitude,1),'YTicklabel',StimResponseFiring.frequencies)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')


    % modulation index
    subplot(2,2,2);
    imagesc(modulation_index(:,:, cluster),'AlphaData', NaNindex_modulation(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = 'multisensory modulation';

    set(gca, 'YDir', 'normal', 'XTick',1:size(modulation_index,2),'XTickLabel',umN,'YTick',1:size(modulation_index,1),'YTicklabel',StimResponseFiring.frequencies)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')

    % additivity index
    subplot(2,2,3);
    imagesc(additivity_index(:,:, cluster),'AlphaData', NaNindex_additivity(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = 'additivity index';

    set(gca, 'YDir', 'normal', 'XTick',1:size(additivity_index,2),'XTickLabel',umN,'YTick',1:size(additivity_index,1),'YTicklabel',StimResponseFiring.frequencies)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')

    % preference index
    subplot(2,2,4);
    imagesc(preference_index(:,:, cluster),'AlphaData', NaNindex_preference(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = 'preference index';

    set(gca, 'YDir', 'normal', 'XTick',1:size(preference_index,2),'XTickLabel',umN,'YTick',1:size(preference_index,1),'YTicklabel',StimResponseFiring.frequencies)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')

end


%% preference index
% % select all units responding to sound and vibration
% unitResponses = unitResponses_all;
% sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
% %sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
% all_stim = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
% %all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);
% %non = table2array(~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SA & ~unitResponses.SO);
% non = table2array(~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SA & ~unitResponses.SO & ~unitResponses.SOM);
% 
% % select and format data
% %sound_vibro_idx = max(all_stim, sound_vibrotac);
% %data = data_all(:,:,:,sound_vibro_idx);
% %aud_data = abs(squeeze(data(1,1,2,:)));
% %som_data = abs(squeeze(mean(data(2:7,3,3,:), 1)));
% %data = data_all(:,:,:,~non);
% %data = squeeze(mean(mean(data, 1, "omitnan"), 2, "omitnan"));
% %data(5,:) = abs(data(2,:)) + abs(data(3,:));
% 
% % pref_index = (dFRsom - dFRaud) / (dFRsom + dFRaud)
% pref_index = (data(2,:) - data(3,:)) ./ (data(2,:) + data(3,:));
% pref_index = (data(2,:) - data(3,:)) ./ (abs(data(2,:)) + abs(data(3,:)));
% 
% pref_index(isnan(pref_index)) = [];
% mean(pref_index)
% std(pref_index)
% 
% [N, edges] = histcounts(pref_index, -6:0.5:2);
% %[N, edges] = histcounts(modulation_index);
% figure;
% for i = 1:(length(edges)-1)
%     binCenters(:,i) = mean(edges(:,i:i+1));
%     i = i+1;
% end
% 
% bar(binCenters, N, 1, 'FaceColor', [0.5 0.5 0.5])
% xline(0, 'k')
% xlabel('Modality preference')
% ylabel('# units')
% set(gca,'fontsize',18)
% clear binCenters
% %%
% 
% % modulation_index = dFRmulti - (dFRsound + dFRsom) 
% modulation_index = data(4,:) - (data(2,:) + data(3,:));
% modulation_index(isnan(modulation_index)) = 0;
% 
% mean(modulation_index)
% std(modulation_index)
% 
% % plot(sort(modulation_index))
% %yline(0)
% %xlabel('unit');
% %ylabel('modualtion index')
% 
% [N, edges] = histcounts(modulation_index, -3:0.5:6);
% %[N, edges] = histcounts(modulation_index);
% figure;
% for j = 1:(length(edges)-1)
%     binCenters(:,j) = mean(edges(:,j:j+1));
%     j = j+1;
% end
% 
% bar(binCenters, N, 1, 'FaceColor', [0.5 0.5 0.5])
% xline(0, 'k')
% xlabel('Modulation strength')
% ylabel('# units')
% set(gca,'fontsize',18)
% 
% clear binCenters

%% compare firing rate unimodal to multimodal stimulus presentation
% VERSION 1
% makes scatter plot of firing rates per responsive unit between two conditions
% data format:
%   firing_mean(freq, amp, condition, resp_cids)
%   conditions = {"OO", "OA", "SO", "SA"};

% select all units responding to sound and vibration
unitResponses = unitResponses_all;
sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
%sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
all_stim = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
%all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);

% select and format data
sound_vibro_idx = max(all_stim, sound_vibrotac);
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

%% compare firing rate unimodal to multimodal stimulus presentation
% VERSION 2
% makes scatter plot of firing rates per responsive unit between two conditions
% data format:
%   firing_mean(freq, amp, condition, :)
%   conditions = {"OO", "OA", "SO", "SA"};

data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids

% compare multi to sound

uSOM = StimResponseFiring.amplitudes;
uAUD = StimResponseFiring.frequencies;

figure
hold on
pos = 1;
for aud = 2:length(uAUD)
    for som = 2:length(uSOM)

        subplot((length(uAUD)-1),(length(uSOM)-1), pos);
        scatter(squeeze(data(aud,1,:)), squeeze(data(aud,som,:)), 'filled', 'k')
        hold on
        xlim(gca, [-10 100])
        ylim(gca, [-10 100])

        xlabel(['\Delta FR sound (' num2str(uAUD(aud)) 'dbSPL)'])
        ylabel(['\Delta FR (' num2str(uAUD(aud)) 'dbSPL *' num2str(uSOM(som)) 'V)'])

        % add 45 degree angle line
        plot([-10, 100],[ -10, 100], 'k')

        pos = pos+1;
    end
end

% compare multi to som pressure

figure
hold on
pos = 1;
for aud = 2:length(uAUD)
    for som = 2:length(uSOM)

        subplot((length(uAUD)-1),(length(uSOM)-1), pos);
        scatter(squeeze(data(1,som,:)), squeeze(data(aud,som,:)), 'filled', 'k')
        hold on
        xlim(gca, [-10 100])
        ylim(gca, [-10 100])

        xlabel(['\Delta FR pressure (' num2str(uSOM(som)) 'V)'])
        ylabel(['\Delta FR (' num2str(uAUD(aud)) 'dbSPL *' num2str(uSOM(som)) 'V)'])

        % add 45 degree angle line
        plot([-10, 100],[ -10, 100], 'k')

        pos = pos+1;
    end
end

%%  compare firing rate unimodal to unimodal stimuli presentation
% for all units plot dFR for each som and aud stimulus combination
% does not give useful plot

uSOM = StimResponseFiring.amplitudes;
uAUD = StimResponseFiring.frequencies;
data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids

figure
hold on
pos = 1;
for aud = 2:length(uAUD)
    for som = 2:length(uSOM)

        subplot((length(uAUD)-1),(length(uSOM)-1), pos);
        scatter(squeeze(data(aud,1,:)), squeeze(data(1,som,:)), 'filled', 'k')
        hold on
        xlim(gca, [-10 100])
        ylim(gca, [-10 100])
        xlabel([num2str(uSOM(som)) ' (V)']);
        ylabel([num2str(uAUD(aud)) ' dB SPL']);

        %axis([-10 90 -10 90])
        plot([-10, 100],[ -10, 100], 'k')

        pos = pos+1;
    end
end


%% linegraph: plot bbn for multiple animals
uparamA = [0 15 30 45 60];

%data_all = squeeze(data_all);
%data = mean(data_all,2);
%errors = std(data_all,[],2) / sqrt(size(data_all, 2));

% select data in new format
%cluster_index = StimResponseFiring_all.unitResponses{:,5} == 1;
%data_all = squeeze(max(StimResponseFiring_all.firing_mean(:,1,2,cluster_index'), StimResponseFiring_all.firing_mean(:,1,1,cluster_index')));
data_all = squeeze(max(StimResponseFiring.firing_mean(:,1,2,:), StimResponseFiring.firing_mean(:,1,1,:)));
data = mean(data_all,2);

figure;
hold on
plot(uparamA, data_all, 'Color', "#CE87AA")
plot(uparamA, data, 'Color', "#9E0E56", 'LineWidth', 3)

%errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(uparamA)

ylabel('\Delta Firing rate (Hz)')
xlabel('Broadband noise intensity (dB SPL)')
ylim([-50 120])
set(gca,'fontsize',18)

%% linegraph: pressure tuning curve
% work with 4D data format? 
%StimResponseFiring = StimResponseFiring_all;

% variables
uparamA = unique(StimResponseFiring.amplitudes)';
%nAmp = length(uAmp);

umN = round((uparamA * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

% format data
%data = squeeze(data_all)';
data = squeeze(StimResponseFiring.firing_mean(1,:,3,:))';
data(:,1) = squeeze(StimResponseFiring.firing_mean(1,1,1,:))';
index = StimResponseFiring.unitResponses.SO;

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

%% scatter plot: pressure tuning curve

data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids
dataToPlot = squeeze(data(1,:,:));
x = StimResponseFiring.amplitudes';
figure
hold on
%scatter(StimResponseFiring.amplitudes', dataToPlot)
for i = 1:length(x)
    swarmchart(x(i), dataToPlot(:,i), 'filled');

    %swarmchart(x(i), dataToPlot(:,i), 20, 'k', 'filled','XJitterWidth',0.8);
    % swarmchart(x, dataToPlot, 20, 'k', 'filled','XJitterWidth',0.8);
end
xlim(gca, [-0.01 0.34])

figure
%scatter(x, dataToPlot)
plot(x,dataToPlot)
hold on
plot(x, mean(dataToPlot'), 'k', 'LineWidth', 1.5)


%% linegraph: vibrotac tuning curves
% of only responsive units
% work with 4D data format?

%uAmp = unique(stimuli_parameters.Stm.Amplitude);
uparamA = [0;0.1;0.3];
%nAmp = length(uAmp);

umN = round((uparamA * 0.158)*1000); %0.158N/V callibation
%nmN = length(umN);

%uFreq = unique(stimuli_parameters.Stm.SomFreq);
uparamB = [0; 10;20;50;100;200;300;400];
%nFreq = length(uFreq);

% select correct condition from SxA session
condition = 3; % ["OO", "OA", "SO", "SA"]
data_all = StimResponseFiring_all.firing_mean;

%responsive units
units_idx = StimResponseFiring_all.unitResponses.SO;
data = repmat(data_all(1,1,1,units_idx), 1, 3);
data(2:8, 2:3,:, :) = data_all(2:8, 2:3, condition, units_idx);
FRvibrotac_mean = mean(data, 4, "omitnan");
FRvibrotac_se = std(data, [],4) / sqrt(size(data, 4));
FRvibrotac_std = std(data, [],4);

%non responsive units
no_resp_idx = ~StimResponseFiring_all.unitResponses.SO;
data_c = repmat(data_all(1,1,1,no_resp_idx), 1, 3);
data_c(2:8, 2:3,:, :)= data_all(2:8, 2:3, condition, no_resp_idx);
FRvibrotac_mean_c = mean(data_c, 4, "omitnan");
FRvibrotac_se_c = std(data_c, [],4) / sqrt(size(data_c, 4));
FRvibrotac_std_c = std(data_c, [],4);

% plot
figure;
hold on

%errorbar(uFreq, FRvibrotac_mean, FRvibrotac_se, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
for cluster = 1:size(data,4)
    %plot(uFreq , data(:, 1, cluster), 'Color', [0.5, 0.5, 0.5])
    plot(uparamB , data(:, 2, cluster), 'Color', "#8BC9E8")

end

for cluster_c = 1:size(data_c,4)
    plot(uparamB , data_c(:, 2, cluster_c), 'Color', [0.5, 0.5, 0.5])
end

plot(uparamB, FRvibrotac_mean(:, 2), 'LineWidth', 3, 'Color', "#3FA2D4")
%errorbar(uFreq, FRvibrotac_mean, FRvibrotac_std, 'b', 'linestyle', 'none', 'LineWidth', 0.5);

plot(uparamB, FRvibrotac_mean_c(:, 2), 'k', 'LineWidth', 3)
%errorbar(uFreq, FRvibrotac_mean_c, FRvibrotac_std_c, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% format figure
%plot(FRvibrotac_mean(2:nFreq, 2:nAmp));

ax = gca;
xlim = ([0 400]);
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Vibrotactile frequency (Hz)')
%legend(num2str(umN(2:end)))
set(ax,'fontsize',18)
xticks([0 10 20 50 100 200 300 400])
ax.XScale = 'log';


%% heatmap: SxA pressure

% format data
data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids
amps = StimResponseFiring.amplitudes(:,1);
force = round((amps * 0.158)*1000); %0.158N/V callibation
freqs = StimResponseFiring.frequencies(:,1);
cids = StimResponseFiring.unitResponses.Cluster';
mouse = StimResponseFiring.unitResponses.MouseNum';
%OutPath = 'D:\DATA\Processed\temp';
savefile = 0;

for cluster = 1:length(cids)
    figure
    imagesc(data(:,:, cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = '	\Delta spike rate (spikes/sec)';

    set(gca, 'YDir', 'normal', 'XTick',1:size(data,2),'XTickLabel',force','YTick',1:size(data,1),'YTicklabel',freqs')
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    %sgtitle(['unit ' num2str(cids(cluster))])

    if savefile
        % save
        figname = sprintf('M%.2i_cluster%i_PxApost', mouse(cluster), cids(cluster));
        saveas(gcf, fullfile(OutPath, figname));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
    end

end
%% ISI histogram

% amp = 1;
% freq = (uFreq == 0);
% con = 1;

% select unit to plot
%cluster = 2;

%StimResponseFiring = StimResponseFiring_all;

OutPath = 'D:\DATA\Processed\temp';

conditions = StimResponseFiring.conditions;
uparamA = StimResponseFiring.amplitudes;
uparamB = StimResponseFiring.frequencies;
nparamA = length(uparamA);
nparamB = length(uparamB);

for cluster = 12 %1:length(StimResponseFiring.cids)
for con = 1:length(conditions)
    for amp = 1:nparamA

        figure;
        hold on
        pos = 1;

        for freq = 1:nparamB

            % select unique stimulus combination
            if strcmp(stimuli_parameters.Par.Rec, "SxA")
                index = stimuli_parameters.Stm.Amplitude == uparamA(amp) & stimuli_parameters.Stm.SomFreq == uparamB(freq) & strcmp(stimuli_parameters.Stm.MMType, conditions{con});
                figTitle = [num2str(uparamB(freq)) 'Hz ' num2str(uparamA(amp)) 'V'];
                figSGTitle = ['Unit: ' num2str(cids(cluster)) ', condition: ' conditions{con} '(' num2str(uparamA(amp)) 'V)'];
            elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
                index = stimuli_parameters.Stm.Intensity == uparamA(amp) & stimuli_parameters.Stm.Mf == uparamB(freq);
                figTitle = [num2str(uparamB(freq)) 'Hz '];
                figSGTitle = ['Unit: ' num2str(cids(cluster)) ', ' num2str(uparamA(amp)) 'dbSPL'];
            end

            idx_rows = find(index);
            spikeISI = [];
            for i = 1:sum(index)
                tISIspikes = aligned_spikes_ISI{idx_rows(i), cluster};
                spikeISI = [spikeISI; tISIspikes];
            end

            [N, edges] = histcounts(spikeISI, -7:0.05:0);
            subplot(3, 3, pos)
            for i = 1:(length(edges)-1)
                binCenters(:,i) = mean(edges(:,i:i+1));
                % binCenters(:,i) = sqrt(edges(:,i)* edges(:,i+1)); % geometric mean
            end

            bar(binCenters, N, 1, 'FaceColor', 'k')
            xline(log(1/uparamB(freq)),'r--');
            %xline(log(1/uFreq(freq)/2),'b--');
            %xline(log(0.75*1/uFreq(freq)),'k--');
            %xline(log(0.25*1/uFreq(freq)),'g--');
            % stairs(edges(1:end-1), N, 'Color', [0.5 0.5 0.5])
            xticks([log(0.001), log(0.01), log(0.1), log(1)])
            set(gca,'TickDir','out')
            xLab = [0.001, 0.01, 0.1, 1];
            xticklabels(xLab)
            xlabel('ISI (sec)')
            ylabel('# spikes')
            %ylim([0,max(N)+10])
            ylim([0,25])
            %legend;
            title(figTitle)
            clear binCenters

            pos = pos+1;

        end

        sgtitle(figSGTitle)

    end

end

% figname = sprintf('M%.2i_cluster%i_VxApre', animal, cids(cluster));
% saveas(gcf, fullfile(OutPath, figname));
% saveas(gcf, fullfile(OutPath, [figname '.jpg']));


end

%% linegraph: vibrotac tuning curve single units
% of only responsive units
% work with 4D data format?

%uAmp = unique(stimuli_parameters.Stm.Amplitude);
uparamA = [0;0.1;0.3];
nparamA = length(uparamA);
umN = round((uparamA * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

%uFreq = unique(stimuli_parameters.Stm.SomFreq);
uparamB = [10;20;50;100;200;300;400];
nparamB = length(uparamB);

% select correct condition from SxA session
condition = 3; % ["OO", "OA", "SO", "SA"]

%responsive units
%units_idx = unitResponses_all.SO;
unit_idx = 12;
data = StimResponseFiring.firing_mean(2:8, 2:3, condition, unit_idx);

% plot
figure;
hold on

plot(uparamB , data(:, 1), 'LineWidth', 3, 'Color', "#8BC9E8")
plot(uparamB , data(:, 2), 'LineWidth', 3, 'Color', "#3FA2D4")

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
uparamA = unique(StimResponseFiring.amplitudes)';
nparamA = length(uparamA);
umN = round((uparamA * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

%% format data
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
FRpressure_mean(1,1) = StimResponseFiring.firing_mean(1,1,1,12);
FRpressure_mean(1,2:9) = StimResponseFiring.firing_mean(1,2:9,3,12);
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

%% linegraph: noise intensity level tuning curve
uparamA = unique(Stm.Intensity);
nparamA = length(uparamA);
cids = cpos.id;
nClusters = length(cids);

% pressure tuning curve mean of only responsive units
for amp = 1:nparamA
    index = Stm.Intensity == uparamA(amp);
    dfiring_mean(amp, :) = mean(stimulusRate(index, resp_cids));
    dfiring_se(amp,:) = std(stimulusRate(index, resp_cids)) / sqrt(length(resp_cids));
end

figure;
hold on
data = mean(dfiring_mean, 2); %mean(firing_mean, 2);
errors = mean(dfiring_se, 2); %mean(firing_sd, 2);
plot(data, 'LineWidth', 1.25)
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(1:nparamA)
xticklabels(uparamA)

ylabel('df Firing rate (spikes/s)')
xlabel('Broadband noise intensity (dB SPL)')

data = mean(data_all);
errors = std(data_all) / sqrt(size(data_all, 1));
figure;
hold on
plot(data, 'LineWidth', 1.25)
errorbar(1:length(data), data, errors, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
xticks(1:length(uparamA))
xticklabels(uparamA)

ylabel('df Firing rate (spikes/s)')
xlabel('Broadband noise intensity (dB SPL)')
set(gca,'fontsize',14)

%% scatter plot: FSL plotting

%responsive units
sound_resp_units = StimResponseFiring.unitResponses.OA;
vibrotac_resp_units = StimResponseFiring.unitResponses.SO;
multi_resp_units = (StimResponseFiring.unitResponses.SA) | (StimResponseFiring.unitResponses.OA & StimResponseFiring.unitResponses.SO);
index = max(max(sound_resp_units, vibrotac_resp_units), multi_resp_units);

%index = ones(length(cids),1); % all untis
figure; hold on;
x_axis = 1:length(cids(index'));

%control = squeeze(MedFSL(1,1,1,:))*1000;
%scatter(x_axis, control, 'k');
sound = squeeze(MedFSL(1,1,2,index))*1000;
scatter(x_axis, sound, 'r');

vibrotac = squeeze(MedFSL(2:nparamB, 2:nparamA, 3, index))*1000;
vibrotac_med = squeeze(median(vibrotac, 1));
scatter(x_axis, squeeze(vibrotac_med(1, :)), 'c')
scatter(x_axis, squeeze(vibrotac_med(2, :)), 'b')

xlabel('cluster')
ylabel('median FSL (ms)')
xticks(1:cids(index'))
xticklabels(cids(index'))
set(gca, 'YScale', 'log')
ylim([0 20^3])
%legend
legend('sound only', 'med 0.1V', 'med 0.3V');
hold off

% mean FSL for responsive units
% sound = squeeze(MedFSL(1,1,2,index))*1000;
median(sound)
median(vibrotac_med(2,:))

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

% select which session to analyse
session = 4;

% get start + end of session
idx = sessions_TTLs(:,1) == session;
session_time = sessions_TTLs(idx, 3);

% get first and last sample from session
session_start = analog_samples(analog_samples(:, 1) == session_time(1));
session_end = analog_samples(analog_samples(:, 1) == session_time(end));

figure;
plot(analog_samples(session_start:session_end), analog_trace(session_start:session_end)); % analog trace during session 4


%% plot FRA CF and BF over depth
% needs FRA quantification step + cpos
for cluster = 1:length(cids)
    scatter(CF(cluster), cpos.depth(cluster), 'ob')
    hold on
    scatter(BestFreq(cluster), cpos.depth(cluster), 'or')
    if ~isnan(CF(cluster)) && ~isnan(BestFreq(cluster))
        text(7, cpos.depth(cluster), num2str(cpos.id(cluster)))
    end
end
ylabel('relative depth (um)')
xlabel('Frequency (kHz)')
ylim([0 800])
legend('CF', 'BF')

%% plot FRA CF +

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

%% -------------------------- Motor signal  -------------------------- %%
% plot analog motor signal per condition
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
