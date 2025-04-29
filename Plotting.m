%% ---------------------- Load data  ---------------------- %%
% aligned spike data

close all
clearvars

KSPath = 'D:\DATA\EphysRecordingsSorted\M20\ICX\'; % kilosort ephys data
OutPath = 'D:\DATA\Processed\M20\ICX'; % output directory
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M20\'; % stimuli parameters

session = 3;
[cids, stimuli_parameters, aligned_spikes, ~, ~, sessions_TTLs, onsetDelay, ~, clusterinfo] = loadData(OutPath, session, BehaviorPath);

%% add spacing where needed

idx = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
for cluster = 1:length(cids)
    for ii = idx'
        aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} + 0.25;
    end
end

%% plot raster + PSTH of single SOM stimulus type - keep
%saveplots = 0; %0 don't save, 1 save plots in OutPath
% currently plots raster + psth of each vibration freq seperatly

close all
%OutPath = 'D:\DATA\Processed\M10'; % output directory
condition = 'SO';
fig = SOMplotting(stimuli_parameters, aligned_spikes, cids, OutPath, 0, condition);

%% cluster selection new
%stimOrder = readtable('D:\DATA\Processed\M12-16_stimOrder.csv');
%load('D:\DATA\Processed\M10-11-19-20\M10-11-19-20_FRA_pre_UnitResponses.mat');
%OutPath = 'D:\DATA\Processed\M10-11-19-20\figures';

cids = StimResponseFiring_all.unitResponses.Cluster;

% selection non sound vibrotac responsive units
%PxA_SOM = [19 32];% [10 17]; % [19 32]; % selected cids
% selected responsive units
vibrotac_nosound = [10121 10334 10400 10441 10457 11212 11247 11259 19153 19265 19287 20277 20290 20303 20306]'; % vibrotac responsive, during plug
pressure_nosound = [10400 10441 10457 11212 11257 19265 19287 19296 20277 20303 20306]'; % pressure responsive, during plug
resp_cids = unique([vibrotac_nosound; pressure_nosound]); % vibrotac and/or pressure responsive, exclusively with earplug
resp_cids_idx = ismember(cids, resp_cids);

% selection all vibrotactile responsive units
selected_cids = [10121 10328 10330 10334 10382 10386 10400 10421 10423 10441 10457 11210 11212 11217 11247 11257 11259 11263 1939 1990 19153 19265 19287 19296 2046 20225 20276 20277 20287 20290 20303 20306]'; % vibrotac and/or pressure responsive, not exclusively with earplug
%selected_cids_idx = ismember(cids, selected_cids);

% selection sound vibrotac responsive units
nonresp_cids = selected_cids(~ismember(selected_cids, resp_cids));
nonresp_cids_idx = ismember(cids, nonresp_cids);

PxA_SOM = [10 17]; % [19 32]; % selected cids


%all = [resp_cids_idx, nonresp_cids_idx, selected_cids_idx];

%% ---------------------- Channel map plotting  ---------------------- %%

% matches unit position to channel map
[~, channelpos] = plot_in_channel_map(KSPath, cpos.clusterinfo);

%% ---------------------- PostCuration Plotting ---------------------- %%
% plot aligned spikes from PostCuration_BK output files


%% ---------------- USE TO MAKE RASTER OVERVIEW PLOT ---------------- %%
% plotting single sessions - in use
% raster plot of single units
% saves figures to OutPath

for session = 2 %relevant_sessions(1):relevant_sessions(2)
    % load corresponsing files
    [cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

    plotResponses(stimuli_parameters, aligned_spikes, cids, OutPath);

    %close all
end

%% plot single unit - in use
% raster plot + PSTH

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

index = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600);
SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
StimType = stimuli_parameters.Stm.Var25(index);
Var = [SOM_Hz, StimType];
%Var =  [stimuli_parameters.Stm.Var25, stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude];
raster_yinc = [5,10,10];
raster_color = [0, 0, 0];

for cluster = 1:length(cids)
%fig = figure; % rasterplot 
figure;
fig = subplot(3,1,1:2); % rasterplot
ax = gca;

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

fig = subplot(3,1,3); 
hold on
% select groups for hist
OO = stimuli_parameters.Stm.Var25 == 1;
OA = stimuli_parameters.Stm.Var25 == 4;
SO = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600) & stimuli_parameters.Stm.Var25 == 3;
SA = (stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq ~= 600) & stimuli_parameters.Stm.Var25 == 2;

[N, edges] = histcounts(vertcat(aligned_spikes{OO, cluster}), preT:binsize:postT); % OO
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(OO))/binsize)) % spike/s
%plot(edges(1:end-1), ((N/sum(OO))/binsize), 'Color', '#D95319', 'LineWidth',1.5) % spike/s
plot(edges(1:end-1), ((N/sum(OO))/binsize), 'k','LineWidth',2.5) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{OA, cluster}), preT:binsize:postT); % OA
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(OA))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(OA))/binsize), 'Color', "#c21069",'LineWidth',2.5) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{SO, cluster}), preT:binsize:postT); % SO
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(SO))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(SO))/binsize), 'Color', "#3eb5f0",'LineWidth',2.5) % spike/s

[N, edges] = histcounts(vertcat(aligned_spikes{SA, cluster}), preT:binsize:postT); % SA
%histogram('BinEdges', edges, 'BinCounts', ((N/sum(SA))/binsize)) % spike/s
plot(edges(1:end-1), ((N/sum(SA))/binsize), 'Color', "#9475CB",'LineWidth',2.5) % spike/s

%format axis
xlabel('Time (s)')
ylabel('Spike rate (spikes/s)')
xline(xlinerange) % on/off set
legend('control', 'sound only', 'vibrotactile only', 'multimodal', 'Location', 'northeast')

ax2 = gca;
ax2.FontSize = 16;
xlim(ax2,xrange);
%ylim(ax2, [0 60])

title(cids(cluster))

figname = sprintf('M11_S%.2i_%s_cluster_%i', str2num(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures\', figname));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures\', [figname '.jpg']));

end

%% PSTH heatmap of all units - pressure & vibrotactile

vibrotac_nosound = [10121 10334 10400 10441 10457 11212 11247 11259 19153 19265 19287 20277 20290 20303 20306]'; % vibrotac responsive, during plug
pressure_nosound = [10400 10441 10457 11212 11257 19265 19287 19296 20277 20303 20306]'; % pressure responsive, during plug
resp_cids = unique([vibrotac_nosound; pressure_nosound]); % vibrotac and/or pressure responsive, exclusively with earplug

% select units to analyse
clusters_to_plot = resp_cids;
clusters_to_plot_idx = ismember(StimResponseFiring_all.unitResponses.Cluster, resp_cids);

% general parameters
preT  = -0.2;
binsize = 0.01; %10ms
conditions = {'OO', 'OA', 'SO', 'SA'};
PSTHt = [];
PSTH = [];
stimonset = 0; %find(edges == 0);

% select correct parameters
if contains(StimResponseFiring_all.StimType, 'VxA')
    postT = 1.2;
    amp = 0.3; % [0 0.1 0.3]
    freq = 50;
    stimoffset = 0.5;%find(edges == 0.5);
elseif contains(StimResponseFiring_all.StimType, 'PxA')
    postT = 0.3;
    amp = 0.3160; %[0 0.013 0.032 0.063 0.095 0.127 0.19 0.253 0.316]
    audint = 45;
    stimoffset = 0.1;
end

for mouse = 1:size(StimResponseFiring_all.MouseNum,2)

    clear tPSTH

    % load in correct data file
    mousestr = num2str(StimResponseFiring_all.MouseNum(mouse));
    OutPath = ['D:\DATA\Processed\M' mousestr '\ICX'];
    BehaviorPath = ['D:\DATA\Behavioral Stimuli\M', mousestr];
    session = StimResponseFiring_all.session(mouse);
    [cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

    % calculate PSTH per condition
    for condition = 1:length(conditions)
        if contains(StimResponseFiring_all.StimType, 'VxA')
            if condition == 1 || condition == 2
                dataidx = (stimuli_parameters.Stm.Var25 == condition);
            else
                dataidx = (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.SomFreq == freq) & (stimuli_parameters.Stm.Var25 == condition);
            end
        elseif contains(StimResponseFiring_all.StimType, 'PxA')
            if condition == 1 || condition == 2
                dataidx = (stimuli_parameters.Stm.Var25 == condition);
            elseif condition == 3
                dataidx = (stimuli_parameters.Stm.Var25 == condition) & (stimuli_parameters.Stm.Amplitude == amp);
            else
                dataidx = (stimuli_parameters.Stm.Var25 == condition) & (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.AudIntensity == audint);
            end
        end
        if ~sum(dataidx)
            warning('no data selected')
            continue
        end

        % TO DO: add delay in aligned_spikes

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
        spikecounts = spikecounts - mean(baseline,2); % units x bins

        % normalize spikecounts into PSTH = sum spike counts / (number events * bin size)
        tPSTH(:,:,condition) = spikecounts / (size(aligned_spikes(dataidx, cluster), 1) * binsize); % cluster x bins x conditions

    end

    % add 0 to SOM sessions in condition 4
    if size(tPSTH,3) ~=4
        tPSTH(:,:,4) = zeros(size(tPSTH(:,:,1)));
    end

    % concatinate PSTH over all mice
    PSTHt = cat(1,PSTHt,tPSTH);

end

% select in PSTH units
close all
PSTH(:,:,:) = PSTHt(clusters_to_plot_idx,:,:);

% sort PSTH based on strongest response
if ~exist('I', 'var')
    [~, I] = sort(mean(PSTH(:,bin_center > stimonset & bin_center < stimoffset, 3),2), 1); % avg stimperiod firing during vibrotac only
end

PSTH(:,:,:) = PSTH(I,:,:);

% select correct PSTH part
if contains(StimResponseFiring_all.StimType, 'PxA') && condition ~=3
    toplot = PSTH(PxA_SOM(1):PxA_SOM(2),:,:);
    % clusters_idx = 
else
    toplot = PSTH(:,:,:);
end

fontsize = 16;

% plot sound subtracted PSTH
figure;
subplot(3,1,1); % mulitmodal
plot(bin_center, mean(toplot(:,:,4),1), 'k', 'LineWidth', 1)
title('Condition: SA')
if contains(StimResponseFiring_all.StimType, 'VxA')
    sgtitle([num2str(freq) 'Hz, ' num2str(amp) 'V'])
elseif contains(StimResponseFiring_all.StimType, 'PxA')
    sgtitle([num2str(audint) 'dB, ' num2str(amp) 'V'])
end

subplot(3,1,2); % sound only
plot(bin_center, mean(toplot(:,:,2),1), 'k', 'LineWidth', 1)
title('Condition: OA')
subplot(3,1,3); % avg PSTH of multimodal - sound only trace
plot(bin_center,(mean(toplot(:,:,4),1)- mean(toplot(:,:,2),1)), 'k', 'LineWidth', 1)
title('Residual SA - OA')
ylabel('\Delta Firing rate (spikes/s)')
xlabel('time (sec)')

% save
if savefigs
    figname = sprintf('M10-11-19-20_PSTH_SA-OA_Residual_%s', StimResponseFiring_all.StimType);
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', [figname '.jpg']));
end
close(gcf)

for condition = 2:4 % condition = 4; % {'OO', 'OA', 'SO', 'SA'};

    if contains(StimResponseFiring_all.StimType, 'PxA') && condition ~=3
        toplot = PSTH(PxA_SOM(1):PxA_SOM(2),:,condition);
    else
        toplot = PSTH(:,:,condition);
    end

    % plot individual + mean PSTH traces
    figure;
    subplot(2,1,1);
    plot(bin_center,toplot)
    title('individual traces')
    ylabel('\Delta Firing rate (spikes/s)')
    xlabel('time (sec)')
    
    subplot(2,1,2);
    plot(bin_center,mean(toplot,1), 'k', 'LineWidth', 1)
    title('mean PSTH')
    ylabel('\Delta Firing rate (spikes/s)')
    xlabel('time (sec)')

    if contains(StimResponseFiring_all.StimType, 'VxA')
        sgtitle(['Condition: ' conditions{condition} ' (' num2str(freq) 'Hz, ' num2str(amp) 'V)'])
    elseif contains(StimResponseFiring_all.StimType, 'PxA')
        sgtitle(['Condition: ' conditions{condition} ' (' num2str(audint) 'dB, ' num2str(amp) 'V)'])
    end

    if savefigs
        figname = sprintf('M10-11-19-20_PSTH_%s_%s', conditions{condition}, StimResponseFiring_all.StimType);
        saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', figname));
        saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', [figname '.jpg']));
    end
    close(gcf)

    % plot avg PSTH above heatmap
    figure('position',[-1611,418,880,565]);
    fig = tiledlayout(3,1);
    t1 = nexttile(1);
    plot(bin_center,mean(toplot,1), 'k', 'LineWidth', 1.25)

    %plot heatmap
    t2 = nexttile([2 1]);
    
    % accomodate for different numbers of units in each condition
    if  contains(StimResponseFiring_all.StimType, 'PxA') && condition ~=3
        imagesc(edges, 1:size(toplot,1), toplot) % (1:sum(clusters_to_plot(PxA_SOM(1):PxA_SOM(2))))
        %yticks(1:length(selected_cids(PxA_SOM(1):PxA_SOM(2))))
        %yticklabels(selected_cids(I(PxA_SOM(1):PxA_SOM(2))))
    else
        imagesc(edges, (1:sum(clusters_to_plot_idx)), toplot)
        %yticks(1:1:length(selected_cids))
        %yticklabels(selected_cids(I))
    end

    

    if contains(StimResponseFiring_all.StimType, 'VxA')
        t1.YLim = [-2 100];
        clim([-5 250])
        t2.YTick = [1,size(toplot,1)];
        t2.YTickLabel = {num2str(size(toplot,1)), '1'};
        % if condition == 2 % OA
        %     xlim(t1, [-0.2 1.2])
        %     xlim(t2, [-0.2 1.2])
        % elseif condition == 3 % SO
        %     xlim(t1, [-0.2 0.7])
        %     xlim(t2, [-0.2 0.7])
        % elseif condition == 4 % SA
            xlim(t1, [-0.2 1.2])
            xlim(t2, [-0.2 1.2])
        %end
        sgtitle(['Condition: ' conditions{condition} ' (' num2str(freq) 'Hz, ' num2str(amp) 'V)'])
    elseif contains(StimResponseFiring_all.StimType, 'PxA')
        t1.YLim = [-2 80];
        clim([-5 200])
        t2.YTick = [1,size(toplot,1)];
        t2.YTickLabel = {num2str(size(toplot,1)), '1'};
        xlim(t1, [-0.1 0.2])
        xlim(t2, [-0.1 0.2])
        sgtitle(['Condition: ' conditions{condition} ' (' num2str(audint) 'dB, ' num2str(amp) 'V'])
    end

    %format axis
    ylabel(t1, '\Delta spikes/s')
    xlabel(fig, 'Time (s)', 'fontsize',fontsize);
    set(t1,'fontsize',fontsize)
    ylabel(t2, 'units')
    set(t2,'fontsize',fontsize, 'TickDir','out')

    %nexttile([2 1])
    colormap('hot'); % Choose a color map
    cb = colorbar(gca);
    cb.Label.String = '\Delta spike/sec';
    cb.FontSize = fontsize;

    if savefigs
        figname = sprintf('M10-11-19-20_PSTH_heatmap_%s_%s', conditions{condition},StimResponseFiring_all.StimType);
        saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', figname));
        saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', [figname '.jpg']));
    end
end


%% ------------------- dataQuantification Plotting ------------------- %%
% plot aligned spikes from PostCuration_BK output files

OutPath = 'D:\DATA\Processed\M10';
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M10';
session = 2;
[~, ~, ~, ~, ~, ~, ~, StimResponseFiring] = loadData(OutPath, session, BehaviorPath);

%% -------------------------- PLOTTING START -------------------------- %%
%% combine unit responses
unitResponses = unitResponses_all;
sound = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & ~unitResponses.SO);
tactile = table2array(unitResponses.SO & ~unitResponses.OA & ~unitResponses(:,7));
sound_tactile = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO|(unitResponses.SA & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
all_stim = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);
non = table2array(~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SA & ~unitResponses.SO);

clear data dataS unitResponses_all


%% venn diagram: response ratio groups
%sets = ["sound" "vibrotactile" "pressure"];
sets = ["sound" "tactile"];
total = 64; % %size(data_all,4);
sound_resp = 50;
tactile_sound = 9; % 27 for all tac resp
tactile_resp = 17; % 35 for all tac resp
%tactile_only = 8;

%labels = [round((sum(sound_resp)/total)*100), round((sum(tactile_resp)/total)*100), round((sum(tactile_sound)/total)*100)]; % order: A, B, C, A&B, A&C, B&C and A&B&C
labels = [sound_resp-tactile_sound, tactile_resp-tactile_sound, tactile_sound]; % order: A, B, C, A&B, A&C, B&C and A&B&C

venn(2,'sets',sets,'labels',labels, 'colors', [0.808 0.529 0.666;0.247 0.635 0.831], 'alpha',0.75);

%% barplot: average response strength per condition
% data format:
%   firing_mean(freq, amp, condition, :)

conditions = ["OO", "OA", "SO", "SA"];
unitResponses = StimResponseFiring_all.unitResponses;
data_all = StimResponseFiring_all.firing_mean;

% select all units responding to sound and vibration
%sound_vibrotac = table2array(((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO) |unitResponses.SA);
%sound_vibrotac = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO & ~unitResponses.SOM |(unitResponses.SA & ~unitResponses.SOM & ~unitResponses.OA & ~unitResponses(:,7) & ~unitResponses.SO));
%all_stim = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SO);
%all = table2array((unitResponses.OA == 1 | unitResponses(:,7) == 1) & unitResponses.SOM & unitResponses.SO);

% select and format data
%sound_vibro_idx = max(all_stim, sound_vibrotac);
sound_vibro_idx = ismember(StimResponseFiring_all.unitResponses.Cluster, resp_cids);
%data = (mean(data, 1, "omitnan"));
%data = squeeze(mean(mean(data, 1, "omitnan"), 2, "omitnan"));

% make data matrix to plot
data_OO = squeeze(data_all(1,1,1,sound_vibro_idx));
data_OA = squeeze(data_all(1,1,2,sound_vibro_idx)); % if vibrotac
%data_OA = squeeze(data_all(1,1,2,sound_vibro_idx)); % if pressure TO DO
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

data = [data_OO'; data_OA'; data_SO'; data_SA']; % max vibrotac
data(5,:) = data(2,:) + data(3,:);
%data = log(data);
%errors = std(data(:, :), 0, 2); % / sqrt(size(data, 2));

ymin = -5;
ymax = 70;
%%plot bar graph
figure;
set(gcf,'position',[500,150,900,700])
hold on
fig = bar(mean(data,2)); % avg FR of all responsive units
x = repmat((1:5)',1,length(maxValueIdx));
for i = 1:5
    swarmchart(x(i,:), data(i,:), 20, 'k', 'filled','XJitterWidth',0.4);
end

% bar colors
set(fig, 'FaceColor', 'Flat')
fig.CData = [0 0 0; 0.808 0.529 0.666; 0.247 0.635 0.831; 0.580 0.455 0.651; 0.5 0.5 0.5];

% axis
set(gca,'fontsize',18)
%xlabel('conditions')
xticks(1:5)
xticklabels(["control", "sound", "vibrotactile", "multimodal", "sound + vibrotactile"])
ylabel('\Delta Firing rate (spikes/s)')
ylim([ymin ymax])
title(['Condition: ' StimResponseFiring_all.StimType])

hold off

% line graph version
% data format: condition (oo oa so sa oa+so) x units

% TODO: color code responsive units

figure;
subplot(1,2,1)
plot([data_OO,data_OA,data_SA,data_SO]','-o','MarkerSize', 8, 'MarkerFaceColor','k','Color',[.5,.5,.5], 'LineWidth', 1.25)
ylabel('\Delta Firing rate (spikes/s)')
xticks(1:4)
xticklabels(["control", "sound", "multimodal", "vibrotactile"])
ylim([ymin ymax])

% if savefigs
%     figname = sprintf('M10-11-19-20_stimrespchanges_%s_', StimResponseFiring_all.StimType);
%     saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
%     saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
% end

%figure;
subplot(1,2,2)
plot([data_SA,data_OA+data_SO]','-o','MarkerSize', 8, 'MarkerFaceColor','k','Color',[.5,.5,.5], 'LineWidth', 1.25)
%ylabel('\Delta Firing rate (spikes/s)')
xticks(1:2)
xticklabels(["multimodal", "sound + vibrotactile"])
ylim([ymin ymax])

if savefigs
    figname = sprintf('M10-11-19-20_stimrespchanges_additivity_%s_', StimResponseFiring_all.StimType);
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
end



%% heatmap: calculate multisensory integration indexes

% Only works with PxA 9not VxA) stimuli

% basic formulas:
% preference_index = (dFRsom - dFRaud) / (dFRsom + dFRaud)
% additivity_index = (dFRmulti - dFRsom - dFRaud) / (dFRsom + dFRaud)
% enhancement_magnitude = (dFRmulti - dFRmax) / dFRmax
% modulation_index = dFRmulti - (dFRsound + dFRsom)

% TO DO:
% remove data if FR of both modalities < 1Hz
% add FR heatmap

StimResponseFiring = StimResponseFiring_all;

% add monosensory FR to same matrix
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

            % NaN if FR of both modalities < 1Hz
            if abs(dFRsound) <= 1 && abs(dFRvibrotac) <= 1
                dFRsound = NaN;
                dFRvibrotac = NaN;
                dFRmulti = NaN;
                dFRmax = NaN;
            end

            dFRmono_sign = [signvalue(i,1,c), signvalue(1,j,c)];
            dFRmulti_sign = signvalue(i,j,c);
            %dFRmax = max(abs(data(1,j,c)), abs(data(i,1,c)));
            %dFRmax = max(data(1,j,c),data(i,1,c),"ComparisonMethod","abs");

            % calculate multimodal enhancement index for each stim combination
            %enhancement_magnitude(i,j,c) = (dFRmulti - dFRmax) / dFRmax; % multisensory enhancement index
            enhancement_magnitude(i,j,c) = ((((dFRmulti_sign * dFRmono_sign(max_index)) * abs(dFRmulti)) - abs(dFRmax)) / abs(dFRmax))*100; % multisensory enhancement index = (dFRmulti - dFRmax) / dFRmax
            modulation_index(i,j,c) = dFRmulti - (dFRsound + dFRvibrotac); % multisensory modulation index, sign indicated direction multi
            additivity_index(i,j,c) = ((dFRmulti - dFRvibrotac - dFRsound) / abs(dFRvibrotac + dFRsound))*100;
            preference_index(i,j,c) = (dFRvibrotac - dFRsound) / (abs(dFRvibrotac) + abs(dFRsound)); % - pref sound; + pref vibrotac

        end
    end
end

% add NaN to irrelevant combinations (unisensory except control)
enhancement_magnitude(1,2:end,:) = NaN;
enhancement_magnitude(2:end,1,:) = NaN;

modulation_index(1,2:end,:) = NaN;
modulation_index(2:end,1,:) = NaN;

additivity_index(1,2:end,:) = NaN;
additivity_index(2:end,1,:) = NaN;

preference_index(1,2:end,:) = NaN;
preference_index(2:end,1,:) = NaN;

clear i j c

% Plotting multisensory integration measures

%close all

NaNindex_enhancement = ~isnan(enhancement_magnitude) & ~isinf(enhancement_magnitude);
NaNindex_modulation = ~isnan(modulation_index) & ~isinf(modulation_index);
NaNindex_additivity = ~isnan(additivity_index) & ~isinf(additivity_index);
NaNindex_preference = ~isnan(preference_index) & ~isinf(preference_index);

% variables
uAmp = unique(StimResponseFiring.amplitudes)';
umN = round((uAmp * 0.158)*1000); %0.158N/V callibation
uInt = StimResponseFiring.frequencies(:,end);

fontsize = 14;

for cluster = 1:length(StimResponseFiring.unitResponses.Cluster)
    figure('Position',[100,100,1800,500]);
    colormap('parula'); % Choose a color map

    % delta FR
    subplot(1,3,1);
    imagesc(data(:,:, cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = '	\Delta spike rate (spikes/sec)';
    cb.FontSize = fontsize;
    set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN','YTick',1:length(uInt),'YTicklabel',uInt)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    set(gca,'fontsize',fontsize)
    title('Average firing rate')

    % enhancement index
    subplot(1,3,2);
    imagesc(enhancement_magnitude(:,:, cluster),'AlphaData', NaNindex_enhancement(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    %clim([minV, maxV]);
    cb.Label.String = 'multisensory index (%)'; %cb.Label.String = 'multisensory enhancement index';
    cb.FontSize = fontsize;
    set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN,'YTick',1:length(uInt),'YTicklabel',uInt)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    set(gca,'fontsize',fontsize)
    title('Multisensory modulation')

    % additivity index
    subplot(1,3,3);
    imagesc(additivity_index(:,:, cluster),'AlphaData', NaNindex_additivity(:,:,cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = 'additivity index (%)';
    cb.FontSize = fontsize;
    set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN,'YTick',1:length(uInt),'YTicklabel',uInt)
    xlabel('stimulus intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    set(gca,'fontsize',fontsize)
    title('Multisensory additivity')
 
    % save responsive units only
     if ismember(StimResponseFiring.unitResponses.Cluster(cluster), resp_cids)
         sgtitle(['Responsive unit (' num2str(StimResponseFiring.unitResponses.Cluster(cluster)) ')'])
         figname = sprintf('M10-11-19-20_unit-%i_multisensory-indexes_%s', StimResponseFiring.unitResponses.Cluster(cluster), StimResponseFiring.StimType);
         saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
         saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
     else
         sgtitle(['Non responsive unit (' num2str(StimResponseFiring.unitResponses.Cluster(cluster)) ')'])
         close(gcf)
     end


end

% % preference index
% subplot(1,3,1);
% imagesc(preference_index(:,:, cluster),'AlphaData', NaNindex_preference(:,:,cluster))
% cb = colorbar(gca, 'eastoutside');
% cb.Label.String = 'modality preference';
% clim([-1, 1]);
% cb.Ticks = [-1, 0, 1];
%cb.TickLabels = {'tactile preference' 'no preference' 'sound preference'};
% set(gca, 'YDir', 'normal', 'XTick',1:size(preference_index,2),'XTickLabel',umN,'YTick',1:size(preference_index,1),'YTicklabel',StimResponseFiring.frequencies(:,1))
% xlabel('stimulus intensity (mN)')
% ylabel('bbn intensity (dB SPL)')
% % modulation index
% subplot(2,2,2);
% imagesc(modulation_index(:,:, cluster),'AlphaData', NaNindex_modulation(:,:,cluster))
% cb = colorbar(gca, 'eastoutside');
% cb.Label.String = 'multisensory modulation';
% set(gca, 'YDir', 'normal', 'XTick',1:size(modulation_index,2),'XTickLabel',umN,'YTick',1:size(modulation_index,1),'YTicklabel',StimResponseFiring.frequencies(:,1))
% xlabel('stimulus intensity (mN)')
% ylabel('bbn intensity (dB SPL)')
% %% multisensory index tests
% StimResponseFiring = StimResponseFiring_all;
% data = squeeze(sum(StimResponseFiring.firing_mean, 3, 'omitnan')); %nAud x nSom x cids
% unit_idx = 57;
% data2plot = data(:,:,unit_idx);
% dFRmulti = data2plot;
% dFRsound = repmat(dFRmulti(:,1),[1,9]);
% dFRvibrotac = repmat(dFRmulti(1,:),[5,1]);
% dFRmulti_sign = sign(dFRmulti);
% dFRmax = max(dFRsound,dFRvibrotac);
% dFRmax = max(dFRsound, dFRvibrotac,"ComparisonMethod","abs");
% dFRmono_sign = sign(dFRmax);
% 
% msi = (data2plot - dFRmax)./dFRmax * 100; % msi = ((((dFRmulti_sign * dFRmono_sign(max_index)) * abs(dFRmulti)) - abs(dFRmax)) / abs(dFRmax))*100;
% enhancement_magnitude = ((((dFRmulti_sign .* dFRmono_sign) .* abs(dFRmulti)) - abs(dFRmax)) ./ abs(dFRmax))*100; % multisensory enhancement index = (dFRmulti - dFRmax) / dFRmax
% 
% ax1 = subplot(2,3,1); imagesc(data2plot);ax1.YDir='normal';colorbar
% subplot(2,3,2); plot(data2plot(1,:));title('somato')
% subplot(2,3,3); plot(data2plot(:,1));title('aud')
% ax2 = subplot(2,3,4); imagesc(enhancement_magnitude);ax2.YDir='normal';colorbar

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


%% linegraph: plot broadband noise (bbn) for multiple animals
%data_all = [StimResponseFiring_all.firing_mean(1,1,1,selected_cids_idx); StimResponseFiring_all.firing_mean(2:end,1,2,selected_cids_idx)];
%data = squeeze(data_all(:,:,:,PxA_SOM(1):PxA_SOM(2)));

close all
%uInt = unique(Stm.Intensity);
uInt = [0 15 30 45 60];
nInt = length(uInt);

%format data
data = [StimResponseFiring_all.firing_mean(1,1,1,:); StimResponseFiring_all.firing_mean(2:end,1,2,:)];
data = squeeze(data)'; % cids x naudint
FRnoise_resp_mean = mean(data(resp_cids_idx,:), 'omitnan');
FRnoise_nonresp_mean = mean(data(nonresp_cids_idx,:), 'omitnan');

%plot figure
figure;
hold on

% plot selection non sound vibrotac responsive units
plot(uInt, data(resp_cids_idx,:), 'Color',  "#e06ca5", 'LineWidth', 0.1)
plot(uInt, FRnoise_resp_mean, 'Color',  "#c21069", 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med, FRpressure_iqr, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% plot selection sound vibrotac responsive units
plot(uInt, data(nonresp_cids_idx,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
plot(uInt, FRnoise_nonresp_mean, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3)

% format figure
set(gca,'fontsize',16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Broadband noise intensity (dB SPL)')
xticks(uInt)
xticklabels({'-Inf', '15', '30', '45', '60'});
ylim([-5 140])

figname = sprintf('M10-11-19-20_noise tuning curve_%s_rep-vs-nonrep', StimResponseFiring_all.StimType);
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));


% format data
FRnoise_selected_mean = mean(data(selected_cids_idx,:), 'omitnan');
FRnoise_nonselected_mean = mean(data(~selected_cids_idx,:), 'omitnan');

%plot figure
figure;
hold on

% plot selection non sound vibrotac responsive units
plot(uInt, data(selected_cids_idx,:), 'Color',  "#e06ca5", 'LineWidth', 0.1)
plot(uInt, FRnoise_selected_mean, 'Color',  "#c21069", 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med, FRpressure_iqr, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% plot selection sound vibrotac responsive units
plot(uInt, data(~selected_cids_idx,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
plot(uInt, FRnoise_nonselected_mean, 'k', 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med_c, FRpressure_iqr_c, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% format figure
set(gca,'fontsize',16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Broadband noise intensity (dB SPL)')
xticks(uInt)
xticklabels({'-Inf', '15', '30', '45', '60'});
ylim([-5 140])

figname = sprintf('M10-11-19-20_noise tuning curve_%s_selected-vs-noselected', StimResponseFiring_all.StimType);
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));

%% linegraph: single unit pressure

% variables
uInt = StimResponseFiring_all.amplitudes(:,1)';
umN = round((uInt * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

% format data
data = [StimResponseFiring_all.firing_mean(1,1,1,:), StimResponseFiring_all.firing_mean(1,2:end,3,:)];
data = squeeze(data)'; % cids x namp
index = 47;

%plot figure
figure;
hold on
plot(umN, data(index,:), 'Color',  "#267165", 'LineWidth', 3)

% format figure
xlim = ([0 umN(nmN)]);
ylim([-3 30])
set(gca,'fontsize', 16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Stimulus strength (mN)')

%figname = sprintf('M10-11-19-20_pressure tuning curve_%s_resp-vs-nonresp', StimResponseFiring_all.StimType);
%saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
%saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));

%% linegraph: pressure tuning curve

% variables
uInt = StimResponseFiring_all.amplitudes(:,1)';
umN = round((uInt * 0.158)*1000); %0.158N/V callibation
nmN = length(umN);

% format data
data = [StimResponseFiring_all.firing_mean(1,1,1,:), StimResponseFiring_all.firing_mean(1,2:end,3,:)];
data = squeeze(data)'; % cids x namp
FRpressure_resp_mean = mean(data(resp_cids_idx,:), 'omitnan');
FRpressure_nonresp_mean = mean(data(nonresp_cids_idx,:), 'omitnan');

%% plot figure
figure;
hold on
plot(umN, data(resp_cids_idx,:), 'Color',  "#44AA99", 'LineWidth', 0.1)
plot(umN, FRpressure_resp_mean, 'Color',  "#267165", 'LineWidth', 3)
%plot(umN, data(nonresp_cids_idx,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
%plot(umN, FRpressure_nonresp_mean, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3)

% format figure
%xlim([0 umN(nmN)]);
xticks(umN)
ylim([-3 30])
set(gca,'fontsize', 16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Stimulus strength (mN)')

figname = sprintf('M10-11-19-20_pressure tuning curve_%s_resp-vs-nonresp', StimResponseFiring_all.StimType);
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));

hold off

%% format data
selected_cids_idx = 37;
FRpressure_selected_mean = mean(data(selected_cids_idx,:), 'omitnan');
FRpressure_nonselected_mean = mean(data(~selected_cids_idx,:), 'omitnan');

% plot figure
figure;
hold on
plot(umN, data(selected_cids_idx,:), 'Color',  "#44AA99", 'LineWidth', 3)
plot(umN, FRpressure_selected_mean, 'Color',  "#267165", 'LineWidth', 3)
%plot(umN, data(~selected_cids_idx,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
%plot(umN, FRpressure_nonselected_mean, 'k', 'LineWidth', 3)

% format figure
xlim = ([0 umN(nmN)]);
ylim([-3 30])
set(gca,'fontsize', 16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Stimulus strength (mN)')

figname = sprintf('M10-11-19-20_pressure tuning curve_%s_selected_%i', StimResponseFiring_all.StimType, selected_cids_idx);
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));

%% linegraph: vibrotac tuning curves

% variables
uFreq = [8.9;10;20;50;100;200;300;400]; %uFreq = unique(stimuli_parameters.Stm.SomFreq);

% format data 
data = [StimResponseFiring_all.firing_mean(1,1,1,:); StimResponseFiring_all.firing_mean(2:8,3,3,:)]; % freq (high pressure) x cids
data = squeeze(data)'; % cids x freq (high pressure)
resp_cids_idx = 20;
% different groups to plot
FRvibrotac_mean_resp = mean(data(resp_cids_idx,:), "omitnan");
FRvibrotac_mean_nonresp = mean(data(nonresp_cids_idx,:), "omitnan");

% plot
fig= figure;
fig.Position = [100, 100, 544,584];
hold on

plot(uFreq, data(resp_cids_idx,:), 'Color',  "#8BC9E8", 'LineWidth', 3)
plot(uFreq, FRvibrotac_mean_resp, 'Color',  "#52A3CF", 'LineWidth', 3)

%plot(uFreq, data(nonresp_cids_idx,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
%plot(uFreq, FRvibrotac_mean_nonresp, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3)
%errorbar(uFreq, FRvibrotac_mean_c, FRvibrotac_std_c, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

ax = gca;
ylim([-2 50])
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Vibrotactile frequency (Hz)')
set(ax,'fontsize',18)
xticklabels({'0', '10', '20', '50', '100', '200', '300', '400'})
ax.XScale = 'log';
xticks(uFreq)
xticklabels({'', '10', '20', '50', '100', '200', '300', '400'})
xlim([0 400])


figname = sprintf('M10-11-19-20_vibrotactile tuning curve_%s_resp_%i', StimResponseFiring_all.StimType, cids(20));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));

hold off

%% different groups to plot
FRvibrotac_mean_selec = mean(data(selected_cids_idx,:), "omitnan");
FRvibrotac_mean_nonselec = mean(data(~selected_cids_idx,:),  "omitnan");

% plot
fig= figure;
fig.Position = [100, 100, 1400,550];
hold on

plot(uFreq, data(selected_cids_idx,:), 'Color',  "#8BC9E8", 'LineWidth', 0.1)
plot(uFreq, FRvibrotac_mean_selec, 'Color',  "#3FA2D4", 'LineWidth', 3)

plot(uFreq, data(~selected_cids_idx,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
plot(uFreq, FRvibrotac_mean_nonselec, 'k', 'LineWidth', 3)

ax = gca;
ylim([-2 50])
xlim = ([0 400]);
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Vibrotactile frequency (Hz)')
%legend(num2str(umN(2:end)))
set(ax,'fontsize',18)
xticks([0 10 20 50 100 200 300 400])
%ax.XScale = 'log';

figname = sprintf('M10-11-19-20_vibrotactile tuning curve_%s_selected-vs-nonselected', StimResponseFiring_all.StimType);
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));


%% scatter plot: pressure tuning curve
StimResponseFiring = StimResponseFiring_all;
data = squeeze(sum(StimResponseFiring.firing_mean, 4, 'omitnan')); %nAud x nSom x condition x cids
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

    % if savefile
    %     % save
    %     figname = sprintf('M%.2i_cluster%i_PxApost', mouse(cluster), cids(cluster));
    %     saveas(gcf, fullfile(OutPath, figname));
    %     saveas(gcf, fullfile(OutPath, [figname '.jpg']));
    % end

end
%% ISI histogram

%StimResponseFiring = StimResponseFiring_all;
%OutPath = 'D:\DATA\Processed\temp';

conditions = StimResponseFiring.conditions; % OO OA SO SA
uparamA = StimResponseFiring.amplitudes;
uparamB = StimResponseFiring.frequencies;
nparamA = length(uparamA);
nparamB = length(uparamB);
ISICVs = nan(nparamA, nparamB, length(conditions), length(cids)); % most freq ISI %amp,freq,con,cluster

binsize = 0.05; %50ms
%cluster = 1;
for cluster = 1%close all:length(StimResponseFiring.cids)
    for con = 1:length(conditions)
        for amp = 1:nparamA

            figure;
            hold on
            pos = 1;

            for freq = 1:nparamB         

                %spikeISI 4D cell array (amp, freq, con, cluster)
                toplot = spikeISI{amp, freq, con, cluster};
                [N, edges] = histcounts(toplot, -7:binsize:0); %-7: 0.001 sec ISI
                subplot(3, 3, pos)
                for i = 1:(length(edges)-1)
                    binCenters(:,i) = mean(edges(:,i:i+1));
                    % binCenters(:,i) = sqrt(edges(:,i)* edges(:,i+1)); % geometric mean
                end
              
                % plot ISI histograms
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
                ylim([0,30])
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

% plot CVs
figure

for cluster = 1:length(cids)
    hold on
    scatter(0, squeeze(ISICVs(1,1,1,cluster)), 'r')
    scatter(-10, squeeze(ISICVs(1,1,2,cluster)), 'r')
    scatter(uparamB(2:8),squeeze(ISICVs(3,2:8,4,cluster)), 'c') %SA, 0.3V, all freqs, all units
    scatter(uparamB(2:8),squeeze(ISICVs(3,2:8,3,cluster)), 'k') %SO, 0.3V, all freqs, all units

end
ylabel('CV')
legend('cyan: multimodal, black: tactile')
%set(gca,'Xscale','log')

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


%% ----------- FSL plotting vibrotac x sound ----------- %%

% select responsive units
index = ismember(cids, resp_cids); %max(max(sound_resp_units, vibrotac_resp_units), multi_resp_units);

% Initialize the figure
figure;

% Subplot 1: Scatter plot
subplot(2, 1, 1); % 2 rows, 1 column, 1st subplot
hold on;
x_axis = 1:length(cids(index));
logYTicks = [10 100 1000];

% Select data to plot
control = squeeze(StimResponseFiring_all.FSLmed(1,1,1,index)) * 1000;
sound = squeeze(StimResponseFiring_all.FSLmed(1,1,2,index)) * 1000;
vibrotac = squeeze(StimResponseFiring_all.FSLmed(8,3,3,index)) * 1000; % 0.3V 10Hz
multi = squeeze(StimResponseFiring_all.FSLmed(8,3,4,index)) * 1000; % 0.3V 10Hz

% Make scatter plot
scatter(x_axis, control, 'k');
scatter(x_axis, sound, 'r');
scatter(x_axis, vibrotac, 'b');
scatter(x_axis, multi, 'g');

% Set axis
xlabel('Cluster');
ylabel('Median FSL (ms)');
xticks(1:cids(index));
xticklabels(cids(index));
set(gca, 'YScale', 'log');
yticks(logYTicks);
yticklabels(logYTicks);
legend('Control', 'Sound only', 'Vibrotac 0.3V 400Hz', 'Sound + Vibrotac');
title('Scatter Plot of FSL Data');
hold off;

% Subplot 2: Bar graph
subplot(2, 1, 2); % 2 rows, 1 column, 2nd subplot
hold on;

% Replace Inf with NaN
data = [control'; sound'; vibrotac'; multi'];
data(isinf(data)) = NaN;

% Plot bar graph
bar(mean(data, 2, "omitnan")); % Avg FR of all responsive units

% Add scatter points to bar graph
x = repmat((1:(size(data,1)))', 1, size(data, 2));
for i = 1:size(data,1)
    swarmchart(x(i,:), data(i,:), 20, 'k', 'filled', 'XJitterWidth', 0.4);
end

% Customize bar graph
xlabel('Condition');
ylabel('Mean FSL (ms)');
xticks(1:size(data,1));
xticklabels({'Control', 'Sound', 'Vibrotac', 'Multi'});
set(gca, 'YScale', 'log');
yticks(logYTicks);
yticklabels(logYTicks);
title('Bar Graph of Mean FSL');
hold off;

% Adjust layout (optional)
sgtitle('Combined Scatter Plot and Bar Graph');

%% ----------- FSL plotting pressure ----------- %%

%[cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, onsetDelay, StimResponseFiring, clusterinfo] = loadData(OutPath, session, BehaviorPath);

%responsive units
%sound_resp_units = StimResponseFiring.unitResponses.OA;
%vibrotac_resp_units = StimResponseFiring.unitResponses.SO;
%multi_resp_units = (StimResponseFiring.unitResponses.SA) | (StimResponseFiring.unitResponses.OA & StimResponseFiring.unitResponses.SO);
%resp_units = [276, 277, 290, 303, 306];
index = ismember(cids, resp_cids); %max(max(sound_resp_units, vibrotac_resp_units), multi_resp_units);

% Initialize the figure
figure;

% Subplot 1: Scatter plot
subplot(2, 1, 1); % 2 rows, 1 column, 1st subplot
hold on;
x_axis = 1:length(cids(index));
logYTicks = [10 100 1000];

% Select data to plot
control = squeeze(StimResponseFiring_all.FSLmed(1,1,1,index)) * 1000;
pressure = squeeze(StimResponseFiring_all.FSLmed(1,2:end,3,index)) * 1000;

% Make scatter plot
scatter(x_axis, control, 'k', "filled");
for pressure_int = 1:size(pressure, 1)
    scatter(x_axis, pressure(pressure_int,:));
end

% Set axis
xlabel('Cluster');
ylabel('Median FSL (ms)');
xticks(1:cids(index));
xticklabels(cids(index));
set(gca, 'YScale', 'log');
yticks(logYTicks);
yticklabels(logYTicks);
legend('control')
%legend(num2str(StimResponseFiring_all.amplitudes(:,1)));
title('FSL pressure per unit');
hold off;

% Subplot 2: Bar graph
subplot(2, 1, 2); % 2 rows, 1 column, 2nd subplot
hold on;

% Replace Inf with NaN
data = [control'; pressure];
data(isinf(data)) = NaN;

% Plot bar graph
bar(mean(data, 2, "omitnan"), 'FaceAlpha', 0.6); % Avg FR of all responsive units
bar(median(data, 2, "omitmissing"), 'FaceAlpha', 0.6); % med FR of all responsive units

% Add scatter points to bar graph
x = repmat((1:(size(data,1)))', 1, size(data, 2));
for i = 1:size(data,1)
    swarmchart(x(i,:), data(i,:), 20, 'k', 'filled', 'XJitterWidth', 0.4);
end

% Customize bar graph
xlabel('Pressure (V)');
ylabel('Mean FSL (ms)');
xticks(1:size(data,1));
xticklabels(num2str(StimResponseFiring_all.amplitudes(:,1)));
set(gca, 'YScale', 'log');
yticks(logYTicks);
yticklabels(logYTicks);
legend("mean", "median")
title('Mean FSL over pressure intensities');
hold off;

%% ----------- FSL plotting pressure x sound ----------- %%

%[cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, onsetDelay, StimResponseFiring, clusterinfo] = loadData(OutPath, session, BehaviorPath);

%responsive units
%sound_resp_units = StimResponseFiring.unitResponses.OA;
%vibrotac_resp_units = StimResponseFiring.unitResponses.SO;
%multi_resp_units = (StimResponseFiring.unitResponses.SA) | (StimResponseFiring.unitResponses.OA & StimResponseFiring.unitResponses.SO);
%resp_units = [276, 277, 290, 303, 306];
index = ismember(cids, resp_cids); %max(max(sound_resp_units, vibrotac_resp_units), multi_resp_units);

% Initialize the figure
figure; hold on;

% Subplot 1: Scatter plot
%subplot(2, 1, 1);
x_axis = 1:length(cids(index));
logYTicks = [10 100 1000];

% Replace Inf with NaN
tdata = StimResponseFiring_all.FSLmed;
tdata(isinf(tdata)) = NaN;
% initiate data and fill
data = NaN(size(tdata, 1), size(tdata, 2), sum(index));
%data(1,1,:) = tdata(1,1,1,:);
data(1,1,:) = tdata(1,1,1,index) * 1000;
data(2:end, 1, :) = tdata(2:end,1,2,index) * 1000;
data(1, 2:end, :) = tdata(1,2:end,3,index) * 1000;
data(2:end, 2:end, :) = tdata(2:end,2:end,4,index) * 1000;

% Plot bar graph
%bar(mean(data, 2, "omitnan"), 'FaceAlpha', 0.6); % Avg FR of all responsive units
%bar(median(data, 2, "omitmissing"), 'FaceAlpha', 0.6); % med FR of all responsive units

% variables
uAmp = unique(StimResponseFiring_all.amplitudes)';
umN = round((uAmp * 0.158)*1000); %0.158N/V callibation
uInt = StimResponseFiring_all.frequencies(:,end);

fontsize = 14;

% make heatmap of FSL for each unit
for cluster = 1:length(resp_cids)
    %figure('Position',[100,100,1800,500]);
    figure
    colormap('parula'); % Choose a color map

    % FSL
    imagesc(data(:,:, cluster))
    cb = colorbar(gca, 'eastoutside');
    cb.Label.String = 'FSL (ms)';
    cb.FontSize = fontsize;
    set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN','YTick',1:length(uInt),'YTicklabel',uInt, 'ColorScale','log')
    xlabel('pressure intensity (mN)')
    ylabel('bbn intensity (dB SPL)')
    set(gca,'fontsize',fontsize)
    title(['FSL cluster: ' num2str(resp_cids(cluster))])
end

% make heatmap of median FSL

figure
colormap('parula'); % Choose a color map

% FSL
imagesc(median(data(:,:,:),3,'omitmissing'))
cb = colorbar(gca, 'eastoutside');
cb.Label.String = 'FSL (ms)';
cb.FontSize = fontsize;
set(gca, 'YDir', 'normal', 'XTick',1:length(umN),'XTickLabel',umN','YTick',1:length(uInt),'YTicklabel',uInt)
xlabel('pressure intensity (mN)')
clim([0, 50]);
ylabel('bbn intensity (dB SPL)')
set(gca,'fontsize',fontsize)
title('median FSL')


%% ----------- FSL plotting bbn intensity ----------- %%

%[cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, onsetDelay, StimResponseFiring, clusterinfo] = loadData(OutPath, session, BehaviorPath);

%responsive units
%sound_resp_units = StimResponseFiring.unitResponses.OA;
%vibrotac_resp_units = StimResponseFiring.unitResponses.SO;
%multi_resp_units = (StimResponseFiring.unitResponses.SA) | (StimResponseFiring.unitResponses.OA & StimResponseFiring.unitResponses.SO);
%resp_units = [276, 277, 290, 303, 306];
index = ismember(cids, resp_cids); %max(max(sound_resp_units, vibrotac_resp_units), multi_resp_units);
%index(:) = 1;
% Initialize the figure
figure;

% Subplot 1: Scatter plot
subplot(2, 1, 1); % 2 rows, 1 column, 1st subplot
hold on;
x_axis = 1:length(cids(index));
logYTicks = [10 100 1000];

% Select data to plot
control = squeeze(StimResponseFiring_all.FSLmed(1,1,1,index)) * 1000;
sound = squeeze(StimResponseFiring_all.FSLmed(2:end,1,2,index)) * 1000;

% Make scatter plot
scatter(x_axis, control, 'k', "filled");
for bbn_int = 1:size(sound, 1)
    scatter(x_axis, sound(bbn_int,:));
end

% Set axis
xlabel('Cluster');
ylabel('Median FSL (ms)');
xticks(1:cids(index));
xticklabels(cids(index));
set(gca, 'YScale', 'log');
yticks(logYTicks);
yticklabels(logYTicks);
legend('control')
%legend(num2str(StimResponseFiring_all.amplitudes(:,1)));
title('FSL pressure per unit');
hold off;

% Subplot 2: Bar graph
subplot(2, 1, 2); % 2 rows, 1 column, 2nd subplot
hold on;

% Replace Inf with NaN
data = [control'; sound];
data(isinf(data)) = NaN;

% Plot bar graph
bar(mean(data, 2, "omitnan"), 'FaceAlpha', 0.6); % Avg FR of all responsive units
bar(median(data, 2, "omitmissing"), 'FaceAlpha', 0.6); % med FR of all responsive units

% Add scatter points to bar graph
x = repmat((1:(size(data,1)))', 1, size(data, 2));
for i = 1:size(data,1)
    swarmchart(x(i,:), data(i,:), 20, 'k', 'filled', 'XJitterWidth', 0.4);
end

% Customize bar graph
xlabel('BBN intensity (dB SPL)');
ylabel('Mean FSL (ms)');
xticks(1:size(data,1));
xticklabels(num2str(StimResponseFiring_all.frequencies(:,end)));
set(gca, 'YScale', 'log');
yticks(logYTicks);
yticklabels(logYTicks);
legend("mean", "median")
title('Mean FSL over pressure intensities');
hold off;

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


%% plot FRA CF and BF over depth - for 1 mouse
% needs FRA quantification step done + cpos


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

%% plot FRA CF and BF over depth - for multiple mice
% needs FRA quantification step + cpos
savefigs = 1;
cids = StimResponseFiring_all.unitResponses.Cluster';

% CF & BF
figure;
for cluster = 1:length(cids)
    if ~isnan(BF(cluster)) && ~isnan(CF(cluster))
        scatter(CF(cluster), BF(cluster), 'k', 'filled')
        %text(CF(cluster), BF(cluster), num2str(cids(cluster)))
    end
    hold on
end

ylabel('BF (kHz)')
xlabel('CF (kHz)')
ylim([0 70])
xlim([0 70])

if savefigs
    figname = sprintf('M10-11-19-20_FRA_pre_CF-BF');
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
    close gcf
end
%

% plot CF with BW
excBandWith = squeeze(excBand(:,2,:) - excBand(:,1,:));
figure;
for i = 1:5
    subplot(2,3,i)
    scatter(CF, excBandWith(i,:), 'k', 'filled')

    %title([num2str(deltaInt(i)) 'dB'])
    ylabel(['BW: ' num2str(deltaInt(i)) 'dB'])
    xlabel('CF (kHz)')
    ylim([0 16])
    xlim([0 70])

end

if savefigs
    figname = sprintf('M10-11-19-20_FRA_pre_CF-BW');
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
    close gcf
end

% plot Q = CF / BW
excBandWith = squeeze(excBand(:,2,:) - excBand(:,1,:));
Qval = CF' ./ excBandWith(:,:);
figure;
for i = 1:5
    subplot(2,3,i)
    scatter(CF', Qval(i,:), 'k', 'filled')

    ylabel(['Q' num2str(deltaInt(i)) 'dB'])
    xlabel('CF (kHz)')
    ylim([0 30])
    xlim([0 70])
end

if savefigs
    figname = sprintf('M10-11-19-20_FRA_pre_CF-Q');
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
    close gcf
end

%
figure;
boxplot(Qval', 'Colors', [0 0 0])
% Get the handles of the boxes
h = findobj(gca,'Tag','Box');

% Fill each box with gray color
for j = 1:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'), [0.7 0.7 0.7], 'FaceAlpha', 0.8); % Gray color with 50% transparency
end

% Set outlier color to black
hOutliers = findobj(gca, 'tag', 'Outliers');
for i = 1:length(hOutliers)
    hOutliers(i).MarkerEdgeColor = [0.6 0.6 0.6];
end

% Set median color to black
hMedian = findobj(gca, 'tag', 'Median');
for j = 1:length(hMedian)
    hMedian(j).Color = 'k'; % Set median color to black
    hMedian(j).LineWidth = 2.5;
end

xticklabels([0 10 20 30 40])
xlabel('Bw (dB SPL)')
ylabel('Q value')

if savefigs
    figname = sprintf('M10-11-19-20_FRA_pre_CF-Q_boxplot');
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
    saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));
    close gcf
end

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
