%% ---------------------- Load data  ---------------------- %%
close all
clearvars

% set correct directories
KSPath = '\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Jelle\EphysRecordingsSorted\M23\'; % directory kilosort ephys data
OutPath = '\\store\department\neuw\shared\Aaron Wong\Data\ProcessedData\Jelle\Processed\M23\'; % directory to save to
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M23\'; % directory of stimuli parameter files

%% %% ---------------------- Single unit raster plot  ---------------------- %%
% plotting single sessions
% makes raster plot of all single units
% saves figures to OutPath

% experimental session to analyse
session = 14;

[cids, stimuli_parameters, aligned_spikes, ~, ~, sessions_TTLs, onsetDelay, ~, clusterinfo] = loadData(OutPath, session, BehaviorPath);
disp(['Session loaded: ' stimuli_parameters.Par.Rec ' ' stimuli_parameters.Par.SomatosensoryWaveform])

plotResponses(stimuli_parameters, aligned_spikes, cids, OutPath);

%% --------- Single unit PSTH (1 session, 1 stimulus) --------- %%
% plot PSTH of 1 session of specific stimulus combination

amp = 0.3170; % 0.3170, 0.1270V
freq = 20; % 10 20 50 100 200 300 400 500 600 700 800 900
binsize = 0.02;

session = 14;
[cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

% define shared x-lim parameters
preT  = -0.2;
postT = 0.7;
xrange = [preT, postT];
start_som = 0;
end_som = 0.5;
xlinerange = [start_som end_som];

% make histogram / PSTH
for cluster = 1:length(cids)
        
    figure;
    hold on

    % select groups for histogram
    OO = (stimuli_parameters.Stm.SomFreq == 0) & (stimuli_parameters.Stm.Amplitude == amp);
    SO = (stimuli_parameters.Stm.SomFreq == freq) & (stimuli_parameters.Stm.Amplitude == amp);

    [N, edges] = histcounts(vertcat(aligned_spikes{OO, cluster}), preT:binsize:postT); % OO
    %histogram('BinEdges', edges, 'BinCounts', ((N/sum(OO))/binsize)) % spike/s
    plot(edges(1:end-1), ((N/sum(OO))/binsize), 'k','LineWidth',2.5) % spike/s

    [N, edges] = histcounts(vertcat(aligned_spikes{SO, cluster}), preT:binsize:postT); % SO
    %histogram('BinEdges', edges, 'BinCounts', ((N/sum(SO))/binsize)) % spike/s
    plot(edges(1:end-1), ((N/sum(SO))/binsize), 'Color', "#3eb5f0",'LineWidth',2.5) % spike/s

    %format axis
    xlabel('Time (s)')
    ylabel('Spike rate (spikes/s)')
    xline(xlinerange) % on/off set
    legend('control', 'vibrotactile only', 'Location', 'northeast')

    ax = gca;
    ax.FontSize = 16;
    xlim(ax,xrange);
    %ylim(ax, [0 60]);

    title(cids(cluster))

    %figname = sprintf('M11_S%.2i_%s_cluster_%i', str2num(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
    %saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures\', figname));
    %saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures\', [figname '.jpg']));

end

%% ---------  Single unit PSTH (1 session, all stimuli) --------- %%
% plot PSTH of 1 session of all stimulus combination

session = 14;
[cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

% define shared x-lim parameters
preT  = -0.2;
postT = 0.7;
xrange = [preT, postT];
start_som = 0;
end_som = 0.5;
xlinerange = [start_som end_som];
binsize = 0.02;

all_freqs = unique(stimuli_parameters.Stm.SomFreq);
all_amps = unique(stimuli_parameters.Stm.Amplitude);

% make histogram / PSTH
for cluster = 1:length(cids)
    for amp = 1:length(all_amps) %each amp
        figure;
        hold on
        title(['Cluster: ' num2str(cids(cluster)) ' ' num2str(all_amps(amp)) '(V)'])

        for freq = 1:length(all_freqs) % each freq

            % select groups for histogram
            SO = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == all_amps(amp));
            [N, edges] = histcounts(vertcat(aligned_spikes{SO, cluster}), preT:binsize:postT); % SO
            %histogram('BinEdges', edges, 'BinCounts', ((N/sum(SO))/binsize)) % spike/s

            if all_freqs(freq) == 0 % plot control in black
                plot(edges(1:end-1), ((N/sum(SO))/binsize), 'k','LineWidth',2.5) % spike/s
            else
                % plot(edges(1:end-1), ((N/sum(SO))/binsize), 'Color', "#3eb5f0",'LineWidth',2.5) % spike/s
                plot(edges(1:end-1), ((N/sum(SO))/binsize),'LineWidth',2.5) % spike/s
            end

        end

        %format axis
        xlabel('Time (s)')
        ylabel('Spike rate (spikes/s)')
        legend

        ax = gca;
        ax.FontSize = 16;
        xlim(ax,xrange);
        %ylim(ax, [0 60]);
        xline(xlinerange) % on/off set

    end

end

%% --------- Single unit PSTH (multiple sessions, 1 stimulus) --------- %%
% plot PSTH of several sessions of specific stimulus combination
close all

% select stimuli
amp = 0.3170; % 0.3170, 0.1270
freq = 20;

sessions = [14, 15, 16]; % in order of plotting
cids = loadData(OutPath, sessions(1), BehaviorPath);

% define shared x-lim parameters
preT  = -0.2;
postT = 0.7;
xrange = [preT, postT];
start_som = 0;
end_som = 0.5;
xlinerange = [start_som end_som];
binsize = 0.02;

% make histogram / PSTH
for cluster = 1:length(cids)

    figure;
    hold on

    for session = 1:length(sessions) %each session

        % load data
        [~, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~, ~] = loadData(OutPath, sessions(session), BehaviorPath);

        % select groups for histogram
        OO = (stimuli_parameters.Stm.SomFreq == 0) & (stimuli_parameters.Stm.Amplitude == amp); % control condition
        SO = (stimuli_parameters.Stm.SomFreq == freq) & (stimuli_parameters.Stm.Amplitude == amp); % experimental condition

        [N, edges] = histcounts(vertcat(aligned_spikes{OO, cluster}), preT:binsize:postT); % OO
        %histogram('BinEdges', edges, 'BinCounts', ((N/sum(OO))/binsize)) % spike/s
        plot(edges(1:end-1), ((N/sum(OO))/binsize), 'k','LineWidth',2.5) % spike/s

        [N, edges] = histcounts(vertcat(aligned_spikes{SO, cluster}), preT:binsize:postT); % SO
        %histogram('BinEdges', edges, 'BinCounts', ((N/sum(SO))/binsize)) % spike/s
        %plot(edges(1:end-1), ((N/sum(SO))/binsize), 'Color', "#3eb5f0",'LineWidth',2.5) % spike/s
        plot(edges(1:end-1), ((N/sum(SO))/binsize), 'LineWidth',2.5) % spike/s

    end

    % format figure
    title(['Cluster: ' num2str(cids(cluster))])
    xline(xlinerange) % on/off set
    xlabel('Time (s)')
    ylabel('Spike rate (spikes/s)')
    legend('control location S14', 'location S14', 'control location S15', 'location S15', ...
        'control location S16', 'location S16', 'Location', 'northeast') % EDIT
    %legend

    ax = gca;
    ax.FontSize = 16;
    xlim(ax,xrange);
    %ylim(ax, [0 60]);

end

%% --------- Single unit stimulus + PSTH + raster  --------- %%
% plots raster + psth of each vibration freq in seperate figure
close all

saveplots = 0; %0 don't save, 1 save plots in OutPath
cidtoplot = 260; % choose cluster to plot

condition = 'SO';
fig = SOMplotting(stimuli_parameters, aligned_spikes, cidtoplot, OutPath, saveplots, condition);

%% --------------------- Channel map  --------------------- %%
% matches unit position to channel map
[figa, figb, channelpos] = plot_in_channel_map(KSPath, clusterinfo);


%% --------- Single unit raster - (1 stimulus, different sessions) --------- %%
% plot raster of several sessions of specific stimulus combination

% define shared x-lim parameters
preT  = -0.2;
postT = 0.7;

% select stimuli
amp = 0.3170; % 0.3170, 0.1270
freq = 20;

% select sessions
sessions = [14, 15, 16]; % in order of plotting
cids = loadData(OutPath, sessions(1), BehaviorPath);

SOM_Freq = [];
SOM_Amp = [];
sessions_var = [];
taligned_spikes = [];

% organize variables (Var) and aligned spikes
for session = 1:length(sessions) %each session

    % load data
    [~, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~, ~] = loadData(OutPath, sessions(session), BehaviorPath);

    % select groups for histogram
    %OO = (stimuli_parameters.Stm.SomFreq == 0) & (stimuli_parameters.Stm.Amplitude == amp); % control condition index
    %SO = (stimuli_parameters.Stm.SomFreq == freq) & (stimuli_parameters.Stm.Amplitude == amp); % experimental condition index

    % define stimulus variable space
    %index = (stimuli_parameters.Stm.SomFreq == freq) & (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.AudIntensity == -Inf); % select only high pressure case
    index = (stimuli_parameters.Stm.SomFreq == freq | stimuli_parameters.Stm.SomFreq == 0) & (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.AudIntensity == -Inf); % select only high pressure case

    %linecolor = '#52A3CF';

    tSOM_Freq = stimuli_parameters.Stm.SomFreq(index);
    %tSOM_Amp = stimuli_parameters.Stm.Amplitude(index);
    %Var = [SOM_Hz, SOM_Amp];
    % append Var
    SOM_Freq = [SOM_Freq; tSOM_Freq];
    %SOM_Amp = [SOM_Amp; tSOM_Amp];
    sessions_var = [sessions_var; (repmat(sessions(session), length(tSOM_Freq), 1))];

    % append aligned_spikes
    taligned_spikes = [taligned_spikes; aligned_spikes(index, :)];

end

% rasterplot format: general figure formatting variables
start_stim = 0;
end_stim = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
xlinerange = [start_stim end_stim];
xrange = [preT, postT];

% make raster plot for each unit
for cluster = 1:length(cids)

    % plot both in 1 figure
    figure; hold on
    ax0 = gca;
    Var = [sessions_var, SOM_Freq];

    [f, YTick, ~, ~, ~, YTickLim] = plotraster(gca, taligned_spikes(:, cluster), Var, [0, 0, 0], [10, 10], 1);
    yticks(YTick{2});
    yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
    ylim(f,yrange);
    xlim(f,xrange);

    % add demarcation lines
    xline(xlinerange) % on/off set
    horizontalLine(YTickLim, ax0) % between categories

    % format axis
    xlabel('Time (s)')
    yticklabels({'loc 1', 'loc 1', 'loc 2', 'loc 2', 'loc 3', 'loc 3'})
    ylabel('stimulation location')
    %fig.FontSize = 11;
    title(['Cluster ' num2str(cids(cluster))])

    hold off

    % plot experimental and control conditions in seperate subplots
    figure; hold on
    ax = subplot(2,1,1);
    Var = [SOM_Freq, sessions_var];


    % make experimental condition rasterplot
    [f, YTick, ~, ~, ~, YTickLim] = plotraster(ax, taligned_spikes((SOM_Freq == freq), cluster), Var(SOM_Freq == freq,:), [0, 0, 0], [], 1);

    yticks(YTick{2});
    yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
    ylim(f,yrange);
    xlim(f,xrange);

    % add demarcation lines
    xline(xlinerange) % on/off set
    horizontalLine(YTickLim, ax) % between categories

    % format axis
    %xlabel('Time (s)')
    yticklabels({'loc 1', 'loc 2', 'loc 3'})
    ylabel('stimulation location')
    %fig.FontSize = 11;
    title(['Cluster ' num2str(cids(cluster))])

    % make control condition rasterplot
    ax2 = subplot(2,1,2);
    [f, YTick, ~, ~, ~, YTickLim] = plotraster(ax2, taligned_spikes((SOM_Freq == 0), cluster), Var(SOM_Freq == 0,:), [0, 0, 0], [], 1);

    yticks(YTick{2});
    yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
    ylim(f,yrange);
    xlim(f,xrange);

    % add demarcation lines
    xline(xlinerange) % on/off set
    horizontalLine(YTickLim, ax2) % between categories

    % format axis
    %xlabel('Time (s)')
    yticklabels({'loc 1', 'loc 2', 'loc 3'})
    ylabel('stimulation location')
    xlabel('Time (s)')

    % save plot
    %figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
    %saveas(gcf, fullfile(OutPath, [figname '.jpg']));
    %saveas(fig, fullfile(OutPath, figname));

    hold off
end


%% --------- Population PSTH - (1 stimulus combination, 1 sessions) --------- %%

% DATA & PARAMETER SELECTION
session = 14;
[cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

cidstoplot = cids; % select units of interest

% PSTH parameters
preT  = -0.2;
postT = 0.7;
xrange = [preT, postT];
binsize = 0.02;

freq = 20;
amp = 0.3170;

index = (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.SomFreq == freq);

% ANALYSIS
for cluster = 1:length(cidstoplot)
    [N, edges] = histcounts(vertcat(aligned_spikes{index, cluster}), preT:binsize:postT);

    spikecounts(cluster, :) = N; % cluster x bin
    bin_edges(cluster, :) = edges;
end

% baseline subtraction on spikecounts
bin_center = edges(1:end-1)+0.5*binsize;
baseline = spikecounts(:,bin_center<0); %200ms baseline
spikecounts = spikecounts - mean(baseline,2); % units x bins

% normalize spikecounts into PSTH = sum spike counts / (number events * bin size)
tPSTH(:,:) = spikecounts / (size(aligned_spikes(index, cluster), 1) * binsize); % cluster x bins x conditions

% PLOTTING
figure;
sgtitle([' Stimulation location: ' stimuli_parameters.Par.SomatosensoryLocation])

% individual traces
%subplot(2,1,1)
plot(bin_center,tPSTH)
xlabel('Time (s)')
ylabel('\Delta Firing rate (spikes/s)')
hold on;
plot(bin_center,mean(tPSTH,1), 'k', 'LineWidth', 2.5)
ax1 = gca;
ylim(ax1, [0 60])
ax1.FontSize = 16;

% avg PSTH
%subplot(2,1,1)
%plot(bin_center,mean(tPSTH,1), 'k', 'LineWidth', 1)
%title('mean PSTH')
%xlabel('Time (s)')
%ylabel('\Delta Firing rate (spikes/s)')
%xline(xlinerange) % on/off set
%legend('control', 'sound only', 'vibrotactile only', 'multimodal', 'Location', 'northeast')

%% PSTH heatmap of all units - pressure & vibrotactile
% TO DO: CHANGE FOR JELLE

% DATA SELECTION
resp_cids = [194 180]; % vibrotac and/or pressure responsive, exclusively with earplug
mice_nr = [21 22 23 24 27]; % mice to plot data from
session_nr = [2 2 2 2 2]; % session nr, same order as mice

% PARAMETERS
preT  = -0.2;
postT = 0.7;
binsize = 0.01; %10ms bins
PSTHt = [];
PSTH = [];
stimonset = 0; %find(edges == 0);
stimoffset = 0.5; %find(edges == 0.5);

amp = 0.3; % 0 0.1 0.3
freq = 50; % 10 20 50 100 200 ... 600
all_cids = [];

% ANALYSIS
for mouse = 1:length(mice_nr)

    clear tPSTH

    % load in correct data file
    mousestr = num2str(mice_nr(mouse));
    OutPath = ['D:\DATA\Processed\M' mousestr '\ICX']; % change paths
    BehaviorPath = ['D:\DATA\Behavioral Stimuli\M', mousestr]; % change paths
    session = session_nr(mouse);
    [cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

    % calculate PSTH per condition
    dataidx = (stimuli_parameters.Stm.Amplitude == amp) & (stimuli_parameters.Stm.SomFreq == freq);

    if ~sum(dataidx)
        warning('no data selected')
        continue
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
    spikecounts = spikecounts - mean(baseline,2); % units x bins

    % normalize spikecounts into PSTH = sum spike counts / (number events * bin size)
    tPSTH(:,:) = spikecounts / (size(aligned_spikes(dataidx, cluster), 1) * binsize); % cluster x bins

    % concatinate PSTH over all mice
    PSTHt = cat(1,PSTHt,tPSTH);

    all_cids = [all_cids, cids];

end

close all

% select PSTH units to analyse
clusters_to_plot_idx = ismember(all_cids, resp_cids);
PSTH(:,:) = PSTHt(clusters_to_plot_idx,:);

% sort PSTH based on strongest response
if ~exist('I', 'var')
    [~, I] = sort(mean(PSTH(:,bin_center > stimonset & bin_center < stimoffset, 3),2), 1); % avg stimperiod firing during vibrotac only
end

% select correct PSTH part, for plotting
toplot = PSTH(I,:);

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

sgtitle(['Condition: ' conditions{condition} ' (' num2str(freq) 'Hz, ' num2str(amp) 'V)'])

%     figname = sprintf('M21-22-23-24-27_PSTH_%s_%s', conditions{condition}, StimResponseFiring_all.StimType);
%     saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', figname));
%     saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', [figname '.jpg']));
% close(gcf)

% plot avg PSTH above heatmap
figure('position',[-1611,418,880,565]);
fig = tiledlayout(3,1);
t1 = nexttile(1);
plot(bin_center,mean(toplot,1), 'k', 'LineWidth', 1.25)

%plot heatmap
t2 = nexttile([2 1]);
imagesc(edges, (1:sum(clusters_to_plot_idx)), toplot)

% FORMAT FIGURE
fontsize = 16;
t1.YLim = [-2 100];
clim([-5 250])
t2.YTick = [1,size(toplot,1)];
t2.YTickLabel = {num2str(size(toplot,1)), '1'};
xlim(t1, [-0.2 1.2])
xlim(t2, [-0.2 1.2])
sgtitle(['Condition: ' conditions{condition} ' (' num2str(freq) 'Hz, ' num2str(amp) 'V)'])

ylabel(t1, '\Delta spikes/s')
xlabel(fig, 'Time (s)', 'fontsize',fontsize);
set(t1,'fontsize',fontsize)
ylabel(t2, 'units')
set(t2,'fontsize',fontsize, 'TickDir','out')
colormap('hot'); % Choose a color map
cb = colorbar(gca);
cb.Label.String = '\Delta spike/sec';
cb.FontSize = fontsize;

% figname = sprintf('M21-22-23-24-27_PSTH_heatmap_%s_%s', conditions{condition},StimResponseFiring_all.StimType);
% saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', figname));
% saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20', [figname '.jpg']));

%% linegraph: plot broadband noise (bbn) for multiple animals

close all

% PARAMETERS
uInt = StimResponseFiring_all.amplitudes(:,1); % [0 15 30 45 60];
uInt(isinf(uInt)) = 0; % replace -Inf with 0
nInt = length(uInt);

% DATA SELECTION
resp_cids = [194 180]; % vibrotac and/or pressure responsive, exclusively with earplug
all_cids = StimResponseFiring_all.unitResponses.Cluster; 
resp_cids_idx = ismember(all_cids, resp_cids);

data = squeeze(firing_mean)'; % cids x naudint
FRnoise_resp_mean = mean(data(resp_cids_idx,:), 'omitnan');
%FRnoise_nonresp_mean = mean(data(nonresp_cids_idx,:), 'omitnan');

% PLOTTING
figure;
hold on

% plot selection non sound vibrotac responsive units
plot(uInt, data(resp_cids_idx,:), 'Color',  "#e06ca5", 'LineWidth', 0.1)
plot(uInt, FRnoise_resp_mean, 'Color',  "#c21069", 'LineWidth', 3)
%errorbar(uAmp, FRpressure_med, FRpressure_iqr, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

% plot selection sound vibrotac responsive units
%plot(uInt, data(nonresp_cids_idx,:), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.1)
%plot(uInt, FRnoise_nonresp_mean, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3)

% FORMAT FIGURE
set(gca,'fontsize', 16)
ylabel('\Delta Firing rate (spikes/s)')
xlabel('Broadband noise intensity (dB SPL)')
ylim([-5 85])
 
% figname = sprintf('M10-11-19-20_noise tuning curve_%s', StimResponseFiring_all.StimType);
% saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', figname));
% saveas(gcf, fullfile('D:\DATA\Processed\M10-11-19-20\figures', [figname '.jpg']));

%% ----------------------- LOCAL FUNCTIONS ----------------------- %%
function horizontalLine(YTickLim, fig)

    for i = 1:size(YTickLim,1)
        yline(fig,YTickLim(i,1)-3,':k');
        yline(fig,YTickLim(i,2)+3,':k');
    end

    % % box
    % x = -0.2;w=1;
    % for i = 1:size(YTickLim,1)
    %     r(i) = rectangle(fig,'Position',[x,YTickLim(i,1)-0.5,w,YTickLim(i,2)-YTickLim(i,1)+1],...
    %         'FaceColor',[0,0,0,0],'EdgeColor','none');
    % end
end
