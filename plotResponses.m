function fig = plotResponses(stimuli_parameters, aligned_spikes, cids, OutPath)
% plot data: rater & PSTH

% check input
if isempty(aligned_spikes) || isempty (stimuli_parameters) % to add: look for saved file in folder
    error('no file found')
end

% plot FRA data
if strcmp(stimuli_parameters.Par.Rec, 'FRA')
    for cluster = 1:length(cids)

        fig = figure;
        ax = gca;

        % make rasterplot
        %Var =  [stimuli_parameters.Stm.Intensity,stimuli_parameters.Stm.Freq];
        Var =  [stimuli_parameters.Stm.Freq, stimuli_parameters.Stm.Intensity];
        [f, YTick, YTickLab, varargout] = plotraster(ax, aligned_spikes(:, cluster), Var, [0,0,0], [], 1);
        UFreq = unique(stimuli_parameters.Stm.Freq);
        nFreq = length(UFreq);
        yticks(YTick{1}([1,2:4:nFreq])); yticklabels(UFreq([1,2:4:nFreq]));
        yrange = [min(YTick{end}) - 100, max(YTick{end}) + 100];
        ylim(f,yrange)
        xlim(f,[-.1,.19]);
        set(gcf,'position',[500,150,600,400])

        % format and save
        %title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set])
        sgtitle(['Pure tones raster (unit ' num2str(cids(cluster)) ')' ])
        xlabel('Time (s)')
        ylabel('Stimulus frequency (kHz)')

        % save plot
        figname = sprintf('M%.2i_S%.2i_%s_cluster_%i_raster', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        saveas(fig, fullfile(OutPath, figname));

    end
end

% plot SOM vibrotactile data
if strcmp(stimuli_parameters.Par.Rec, 'SOM')
    %if strcmp(stimuli_parameters.Par.SomatosensoryActuator, 'DC Motor') || strcmp(stimuli_parameters.Par.Rec, 'SOM')

    xrange = [-.2, +.7];

    for cluster = 1:length(cids)

        fig = figure;
        %fig = subplot(2,1,1); % rasterplot
        %set(gcf,'position',[500,150,800,600])

        % define stimulus variable space
        if max(stimuli_parameters.Stm.SomFreq) == 0
            Var = stimuli_parameters.Stm.Amplitude;
        else
            Var = [stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude];
        end

       % [f, YTick, YTickLab] = plotraster(gca, aligned_spikes(:, 1), Var, [0, 0, 0], [], 1);
        % make rasterplot
        [f, YTick, YTickLab, varargout] = plotraster(gca, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
        yticks(YTick{1});
        yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
        ylim(f,yrange);
        yticklabels(unique(stimuli_parameters.Stm.SomFreq));
        xlim(f,xrange);

        % format axis
        xlabel('Time (s)')
        ylabel('Vibrotactile stimulation (Hz)')
        %fig.FontSize = 11;
        title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        % % make histogram / PSTH
        % SOM_idx = (stimuli_parameters.Stm.Amplitude ~= 0); % select stimulus trials
        % ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
        % preT  = -0.2;
        % postT = 0.7;
        % binsize = 0.5;
        % 
        % fig = subplot(2,1,2);
        % 
        % [N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
        % % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
        % plot(edges(1:end-1),((N/sum(SOM_idx == 1))/binsize),'Color', '#D95319','LineWidth',2.5)
        % hold on
        % [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
        % % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')
        % plot(edges(1:end-1),((N/sum(SOM_idx == 0))/binsize),'Color', '#0072BD','LineWidth',2.5)

        % %format axis
        % legend('stimulus', 'control')
        % xlabel('Time (s)')
        % ylabel('Spike rate (Hz)')
        % fig.FontSize = 11;
        % xlim(fig,xrange);
        % 
        % sgtitle(['Cluster ' num2str(cids(cluster))])

        % save plot
        figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        saveas(fig, fullfile(OutPath, figname));

        hold off
    end

end

% plot SxA


% plot noise data
if strcmp(stimuli_parameters.Par.Rec, 'AMn')

    preT  = -0.2;
    postT = 1.5;
    xrange = [preT, postT];

    for cluster = 1:length(cids)

        figure;
        fig = subplot(2,1,1); % rasterplot
        set(gcf,'position',[500,150,900,700])
        %set(gcf,'position',[500,150,600,400])

        % make rasterplot
        Var = [stimuli_parameters.Stm.Intensity, stimuli_parameters.Stm.Mf];
        [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
        yticks(YTick{2});
        yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
        ylim(f,yrange);
        yticklabels(unique(stimuli_parameters.Stm.Mf));
        xlim(f,xrange);

        % format axis
        xlabel('Time (s)')
        ylabel('AM frequency (Hz)')
        ax.FontSize = 11;
        %title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        % set variables for PSTH
        AM_15_idx = stimuli_parameters.Stm.Intensity == 15; % select trials
        AM_30_idx = stimuli_parameters.Stm.Intensity == 30;
        AM_45_idx = stimuli_parameters.Stm.Intensity == 45;
        AM_60_idx = stimuli_parameters.Stm.Intensity == 60;
        AM_0_idx = (stimuli_parameters.Stm.Intensity ~= 15) & (stimuli_parameters.Stm.Intensity ~= 30) & ...
            (stimuli_parameters.Stm.Intensity ~= 45) & (stimuli_parameters.Stm.Intensity ~= 60);

        % set range & figure
        fig = subplot(2,1,2);
        binsize = 0.01;

        % make histogram / PSTH
        [N, edges] = histcounts(vertcat(aligned_spikes{AM_60_idx, cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_60_idx))/binsize), 'Color', '#0072BD', 'LineWidth',1.5) % spike/s
        %plot(edges(1:end-1), ((N/sum(AM_60_idx))/binsize), 'Color', '#4FC5D3', 'LineWidth',1.5) % spike/s
        hold on

        [N,edges] = histcounts(vertcat(aligned_spikes{(AM_45_idx), cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_45_idx))/binsize), 'Color', '#D95319', 'LineWidth',1.5) % spike/s
        %plot(edges(1:end-1), ((N/sum(AM_45_idx))/binsize), 'Color', '#B673C8', 'LineWidth',1.5) % spike/s

        [N,edges] = histcounts(vertcat(aligned_spikes{(AM_30_idx), cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_30_idx))/binsize), 'Color', '#7E2F8E', 'LineWidth',1.5) % spike/s

        [N,edges] = histcounts(vertcat(aligned_spikes{(AM_15_idx), cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_15_idx))/binsize), 'Color', '#77AC30', 'LineWidth',1.5) % spike/s

        [N, edges] = histcounts(vertcat(aligned_spikes{AM_0_idx, cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_0_idx))/binsize), 'Color', 'k', 'LineWidth',1.5) % spike/s

        %format axis
        legend('60dB (SPL)', '45dB (SPL)', '30dB (SPL)', '15dB (SPL)', 'no sound')
        xlabel('Time (s)')
        ylabel('Spike rate (Hz)')
        fig.FontSize = 11;
        xlim(fig,xrange);

        sgtitle(['Broadband noise (unit ' num2str(cids(cluster)) ')'])

        figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        saveas(fig, fullfile(OutPath, figname));

        hold off

    end
end
end