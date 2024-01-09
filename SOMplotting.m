%% lines

%YTickLim: first and last line of each block

for i = 1:size(YTickLim,1)
    yline(ax,YTickLim(i,1)-3,'k');
    yline(ax,YTickLim(i,2)+3,'k');
end

% box
x = -0.2;w=1;
for i = 1:size(YTickLim,1)
    r(i) = rectangle(ax,'Position',[x,YTickLim(i,1)-0.5,w,YTickLim(i,2)-YTickLim(i,1)+1],...
        'FaceColor',[0,0,1,0.1],'EdgeColor','none');
end

%% load correct files
%cids = [90 111 124 126 151 159 169 171 182 188 218 232 256 261 264 265 268 278]; % M4

% 250 trials each
stimuli_parameters_som = load('D:\DATA\Behavioral Stimuli\M4\M4_S04_SOM.mat');
stimuli_parameters_ctrl = load('D:\DATA\Behavioral Stimuli\M4\M4_S05_SOM.mat');

aligned_spikes_som = load('D:\DATA\Processed\M04_S04_SOM_AlignedSpikes.mat');
aligned_spikes_som = aligned_spikes_som.SpkT;
aligned_spikes_ctrl = load('D:\DATA\Processed\M04_S05_SOM_AlignedSpikes.mat');
aligned_spikes_ctrl = aligned_spikes_ctrl.SpkT;

%% format data
stimuli_parameters = vertcat(stimuli_parameters_som.Stm, stimuli_parameters_ctrl.Stm);
aligned_spikes = vertcat(aligned_spikes_som, aligned_spikes_ctrl);

OutPath = 'D:\DATA\Processed\M4';

%% plot SOM data based on Location
xrange = [-.202, +.702];

for cluster = 1:length(cids)

    figname = sprintf('M04_SOM_S04-S05_cluster_%i', cids(cluster));

    figure;

    fig = subplot(2,1,1); % rasterplot
    %set(fig,'position',[500,150,800,600]) %set(gcf,'position',[500,150,600,400])

    %iscontrol = ismember(stimuli_parameters.Location,{'right flank, velcro','right flank, poke'});
    iscontrol = ismember(stimuli_parameters.Location,{'velcro flank L air', 'poke right paw air'});
    Var = [stimuli_parameters.Amplitude, iscontrol];
    %Var =  [stimuli_parameters.Stm.Freq,stimuli_parameters.Stm.Intensity];

    % [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
    [f, YTick, YTickLab,~,~,YTickLim] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
    UAmp = unique(stimuli_parameters.Amplitude);
    nAmp = length(UAmp);
    yticks(YTick{2}); %yticks(YTick{1}([1,2:10:nAmp]));
    %yticklabels({' ' 'Actuator off' ' ' 'Actuator on'}); % yticklabels(UAmp([1,2:10:nAmp])); %yticklabels(unique(stimuli_parameters.Stm.Amplitude));
    yticklabels({' ' '' '' ''});
    yrange = [min(YTick{end}) - 25, max(YTick{end}) + 25]; % yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
    ylim(f, yrange);
    xlim(f, xrange);

    for i = 1:size(YTickLim,1)
        yline(fig,YTickLim(i,1)-3, ':k');
        yline(fig,YTickLim(i,2)+3,':k');
    end

    % format axis
    xlabel('Time (s)')
    %ylabel('Actuator position')
    fig.FontSize = 11;
    %title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

    % make histogram / PSTH
    fig = subplot(2,1,2);
    preT  = -0.2;
    postT = 0.7;
    binsize = 0.03; % bit of smoothing psth

    SOM_on_idx = (stimuli_parameters_som.Stm.Amplitude == 8);
    SOM_off_idx = (stimuli_parameters_som.Stm.Amplitude == 0);
    air_on_idx = (stimuli_parameters_ctrl.Stm.Amplitude == 8);
    air_off_idx = (stimuli_parameters_ctrl.Stm.Amplitude == 0);

    % SOM on aligned_spikes_som
    [N, edges] = histcounts(vertcat(aligned_spikes_som{SOM_on_idx, cluster}), preT:binsize:postT);
    % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
    plot(edges(1:end-1),((N/sum(SOM_on_idx == 1))/binsize),'Color', '#00636C','LineWidth',1.5)
    hold on

    % SOM off aligned_spikes_som
    [N,edges] = histcounts(vertcat(aligned_spikes_som{SOM_off_idx, cluster}), preT:binsize:postT);
    % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')
    plot(edges(1:end-1),((N/sum(SOM_on_idx == 0))/binsize),'Color', '#4FC5D3','LineWidth',1.5)

    % air on aligned_spikes_ctrl
    [N, edges] = histcounts(vertcat(aligned_spikes_ctrl{air_on_idx, cluster}), preT:binsize:postT);
    % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
    plot(edges(1:end-1),((N/sum(air_on_idx == 1))/binsize),'Color', '#501574','LineWidth',1.5)

    % air off
    [N, edges] = histcounts(vertcat(aligned_spikes_ctrl{air_off_idx, cluster}), preT:binsize:postT);
    % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
    plot(edges(1:end-1),((N/sum(air_on_idx == 1))/binsize),'Color', '#B673C8','LineWidth',1.5)

    %format axis
    legend('stimulus on', 'stimulus off', 'control stimulus on', 'control stimulus off')
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    fig.FontSize = 11;
    xlim(fig,xrange);

    %sgtitle(['                    Fur stroke (unit ' num2str(cids(cluster)) ')'])
    sgtitle([stimuli_parameters.Location(1) 'unit' num2str(cids(cluster))])

    % % save plot
    saveas(gcf, fullfile(OutPath, [figname '.jpg']));
    saveas(gcf, fullfile(OutPath, figname));

    % display untill button press
    %waitforbuttonpress
    hold off
    % close

end
%% Plot SOM data of 1 session
aligned_spikes = load('D:\DATA\Processed\M1\aligned_spikes_S07.mat');
aligned_spikes = aligned_spikes.aligned_spikes;
stimuli_parameters = load('D:\DATA\Behavioral Stimuli\M1\M1_S07_SOM.mat');

cluster_info = readtable('D:\DATA\EphysRecordingsSorted\M01\M01_S07\kilosort3\cluster_info.tsv','FileType','text'); % info on clusters
cids = cluster_info.cluster_id(strcmp(cluster_info.group,'good'))';

xrange = [-.2,+.55];

for cluster = 1:length(cids)

    figure;
    fig = subplot(2,1,1); % rasterplot
    set(gcf,'position',[500,150,600,400])

    % PSTH related variables
    SOM_idx = (stimuli_parameters.Stm.Amplitude == 10); % select stimulus trials
    ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
    preT  = -0.2;
    postT = 0.55;
    binsize = 0.03;

    % make rasterplot
    Var = stimuli_parameters.Stm.Amplitude;
    [f, YTick, YTickLab] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
    yticks(YTick{1});
    yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
    ylim(f,yrange);
    yticklabels(unique(stimuli_parameters.Stm.Amplitude));
    xlim(f, xrange);

    % format axis
    xlabel('Time (s)')
    ylabel('Stimulus off / on')
    fig.FontSize = 11;

    % make histogram / PSTH
    fig = subplot(2,1,2);
    [N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
    plot(edges(1:end-1),((N/sum(SOM_idx))/binsize),'Color', '#D95319','LineWidth',1.5)
   % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
    hold on
    [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
    plot(edges(1:end-1),((N/sum(ctrl_idx))/binsize),'Color', '#0072BD','LineWidth',1.5)

    %histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')

    legend('stim on', 'stim off');
    sgtitle(['unit ' num2str(cids(cluster))])

end