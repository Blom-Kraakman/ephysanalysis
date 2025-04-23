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


%% --------- Single unit raster - in progress --------- %%
% plot raster of several sessions of specific stimulus combination
% TO DO: ADD DIFF SESSIONS TOGETHER


% define shared x-lim parameters
preT  = -0.2;
postT = 0.7;
xrange = [preT, postT];


start_stim = 0;
end_stim = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
xlinerange = [start_stim end_stim];

for cluster = 1:length(cids)

    fig = figure;
    ax = gca;

    % for each session
    % define stimulus variable space

    index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.3170) & (stimuli_parameters.Stm.AudIntensity == -Inf); % select only high pressure case
    linecolor = '#52A3CF';

    SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
    SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
    %Var = [SOM_Hz, SOM_Amp];
    Var = [SOM_Hz, sessions];

    % [f, YTick, YTickLab] = plotraster(gca, aligned_spikes(:, 1), Var, [0, 0, 0], [], 1);
    % make rasterplot
    [f, YTick, ~, ~, ~, YTickLim] = plotraster(ax, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
    yticks(YTick{1});
    yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
    ylim(f,yrange);
    xlim(f,xrange);

    % add demarcation lines
    xline(xlinerange) % on/off set
    horizontalLine(YTickLim, ax) % between categories

    % format axis
    xlabel('Time (s)')
    yticklabels(yaxislabels)
    ylabel(yaxistext)
    %fig.FontSize = 11;
    title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

    % save plot
    figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
    saveas(gcf, fullfile(OutPath, [figname '.jpg']));
    saveas(fig, fullfile(OutPath, figname));

    hold off
end
