function fig = plotResponses(stimuli_parameters, aligned_spikes, cids, OutPath)
% plot data: rater & PSTH

% check input
if isempty(aligned_spikes) || isempty (stimuli_parameters)% to add: look for saved file in folder
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

        % figname = sprintf('M02_FRA raster_cluster %i', cids(cluster));
        % saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        % saveas(fig, fullfile(OutPath, figname));

        % display untill button press
        %waitforbuttonpress
       % close

    end
end

% plot SOM data
if strcmp(stimuli_parameters.Par.Rec, 'SOM')

    xrange = [-.2, +.7];

    for cluster = 1:length(cids)

        figure;

        fig = subplot(2,1,1); % rasterplot
        set(gcf,'position',[500,150,800,600])
        %set(gcf,'position',[500,150,600,400])

        % make rasterplot
        Var = stimuli_parameters.Stm.Amplitude;
        %[f, YTick, YTickLab] = plotraster(fig, aligned_spikes(1:40, 1), Var, [0, 0, 0], [], 1);

        [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
        yticks(YTick{1});
        yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
        ylim(f,yrange);
        yticklabels(unique(stimuli_parameters.Stm.Amplitude));
        xlim(f,xrange);

        % format axis
        xlabel('Time (s)')
        ylabel('Stimulus off / on')
        fig.FontSize = 11;
        title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        % make histogram / PSTH
        SOM_idx = (stimuli_parameters.Stm.Amplitude == 10); % select stimulus trials
        ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
        preT  = -0.2;
        postT = 0.7;
        binsize = 0.5;

        fig = subplot(2,1,2);

        [N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
        % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
        plot(edges(1:end-1),((N/sum(SOM_idx == 1))/binsize),'Color', '#D95319','LineWidth',2.5)
        hold on
        [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
        % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')
        plot(edges(1:end-1),((N/sum(SOM_idx == 0))/binsize),'Color', '#0072BD','LineWidth',2.5)

        %format axis
        legend('stimulus', 'control')
        xlabel('Time (s)')
        ylabel('Spike rate (Hz)')
        fig.FontSize = 11;
        xlim(fig,xrange);

        sgtitle(['Cluster ' num2str(cids(cluster))])
        % % save plot
        % %figname = sprintf('PSTH_S06_cluster_%i', cluster);
        % %saveas(gcf, fullfile(OutPath, figname));

        % display untill button press
        %waitforbuttonpress
        hold off
     %   close

    end
end

% plot noise data
if strcmp(stimuli_parameters.Par.Rec, 'AMn')

    xrange = [-.2, +.4];

    for cluster = 1:length(cids)

        figure;
        fig = subplot(2,1,1); % rasterplot
        set(gcf,'position',[500,150,900,700])
        %set(gcf,'position',[500,150,600,400])

        % make rasterplot
        Var = stimuli_parameters.Stm.Intensity;
        [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
        yticks(YTick{1});
        yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
        ylim(f,yrange);
        yticklabels(unique(stimuli_parameters.Stm.Intensity));
        xlim(f,xrange);

        % format axis
        xlabel('Time (s)')
        ylabel('Intensity (dB SPL)')
        ax.FontSize = 11;
        %title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        % set range & figure
        % make histogram / PSTH
        AM_45_idx = stimuli_parameters.Stm.Intensity == 45; % select trials
        AM_60_idx = stimuli_parameters.Stm.Intensity == 60;
        AM_0_idx = (stimuli_parameters.Stm.Intensity ~= 45) & (stimuli_parameters.Stm.Intensity ~= 60);

        fig = subplot(2,1,2);
        preT  = -0.2;
        postT = 0.4;
        binsize = 0.01;

        [N, edges] = histcounts(vertcat(aligned_spikes{AM_60_idx, cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_60_idx))/binsize), 'Color', '#4FC5D3', 'LineWidth',1.5) % spike/s
        hold on
        [N,edges] = histcounts(vertcat(aligned_spikes{(AM_45_idx), cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_45_idx))/binsize), 'Color', '#B673C8', 'LineWidth',1.5) % spike/s

        [N, edges] = histcounts(vertcat(aligned_spikes{AM_0_idx, cluster}), preT:binsize:postT);
        plot(edges(1:end-1), ((N/sum(AM_0_idx))/binsize), 'Color', 'k', 'LineWidth',1.5) % spike/s

        %format axis
        legend('60dB (SPL)', ' 45dB (SPL)', 'no sound')
        xlabel('Time (s)')
        ylabel('Spike rate (Hz)')
        fig.FontSize = 11;
        xlim(fig,xrange);

        sgtitle(['Broadband noise (unit ' num2str(cids(cluster)) ')'])

        figname = sprintf('M02_AMn_cluster %i', cids(cluster));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        saveas(fig, fullfile(OutPath, figname));

        % display untill button press
        %waitforbuttonpress
        hold off
        %close

    end
end
end

