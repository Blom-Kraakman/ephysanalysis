%% plot data: rater & PSTH
% select correct behaviour file
% check for:
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));
%relevant_sessions = [17 21];

stim_counter = 1;

% allign spikes to stimuli (per session)
for file = relevant_sessions(1):relevant_sessions(2)

    % load correct trial and initiate stimuli counter
    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);
    NStim = size(stimuli_parameters.Stm, 1); % nr trials to align per stim session
    last_stim = stim_counter + NStim - 1;

    % plot FRA data
    if strcmp(stimuli_parameters.Par.Rec, 'FRA')
        for cluster = 1:length(cids)

            fig = figure;
            ax = gca;

            % make rasterplot
            %Var =  [stimuli_parameters.Stm.Intensity,stimuli_parameters.Stm.Freq];
            Var =  [stimuli_parameters.Stm.Freq,stimuli_parameters.Stm.Intensity];
            [f, YTick, YTickLab] = plotraster(ax, aligned_spikes(stim_counter:last_stim, cluster), Var, [0,0,0], [], 1);
            UFreq = unique(stimuli_parameters.Stm.Freq);
            nFreq = length(UFreq);
            yticks(YTick{1}([1,2:4:nFreq])); yticklabels(UFreq([1,2:4:nFreq]));
            yrange = [min(YTick{end}) - 100, max(YTick{end}) + 100];
            ylim(f,yrange)
            xlim(f,[-.1,.19]);
            set(gcf,'position',[500,150,600,400])

            % format and save
            title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set])
            xlabel('Time (s)')
            ylabel('Stimulus frequency (kHz)')

            % figname = sprintf('M1_S06_cluster%i', cluster);
            % saveas(fig, fullfile(OutPath, figname));

            disp(stimuli_parameters.Par.Rec)
            disp(NStim)
            disp(stim_counter)
            disp(last_stim)

            % display untill button press
            waitforbuttonpress
            close

        end

    % plot SOM data
    elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')

        xrange = [-.2,+.55];

        for cluster = 1:length(cids)
            % PSTH related variables
            SOM_idx = (stimuli_parameters.Stm.Amplitude == 10); % select stimulus trials
            ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
            preT  = -0.2;
            postT = 0.55;
            binsize = 0.01;

            fig = subplot(2,1,1); % rasterplot
            % set(gcf,'position',[500,150,900,700])
            set(gcf,'position',[500,150,600,400])

            % make rasterplot
            Var = stimuli_parameters.Stm.Amplitude;
            [f, YTick, YTickLab] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
            yticks(YTick{1});
            yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
            ylim(f,yrange);
            yticklabels(unique(stimuli_parameters.Stm.Amplitude));
            xlim(f,xrange);

            % format axis
            xlabel('Time (s)')
            ylabel('Stimulus off / on')
            fig.FontSize = 11;

            % make histogram / PSTH
            fig = subplot(2,1,2);
            [N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
            histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
            hold on
            [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
            histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')

            %format axis
            legend('stimulus', 'control')
            xlabel('Time (s)')
            ylabel('Spike rate (Hz)')
            fig.FontSize = 11;
            xlim(fig,xrange);

            sgtitle(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])
            % % save plot
            % %figname = sprintf('PSTH_S06_cluster_%i', cluster);
            % %saveas(gcf, fullfile(OutPath, figname));

            % display untill button press
            waitforbuttonpress
            close

        end

    end

    stim_counter = stim_counter + NStim; % keep track of how many stimuliu have past

end

%% plot AM
% adapted from SOM
% AMStimTime 200, AMPostTime 400, AMPreTime 100

xrange = [-.25,+.5];

for cluster = 1:length(cids)

    % set range & figure
    preT  = -0.2;
    postT = 0.4;
    binsize = 0.01;
    fig = figure;
    ax = gca;

    % make rasterplot
    Var = stimuli_parameters.Stm.Intensity;
    [f, YTick, YTickLab] = plotraster(ax, aligned_spikes(1:60, cluster), Var, [0, 0, 0], [], 1);
    yticks(YTick{1});
    yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
    ylim(f,yrange);
    yticklabels(unique(stimuli_parameters.Stm.Intensity));
    xlim(f,xrange);

    % format axis
    xlabel('Time (s)')
    ylabel('Intensity dB SPL')
    ax.FontSize = 11;
    title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

    % display untill button press
    waitforbuttonpress
    close

end

%% plot SOM
% raster & PSTH
% to do: select correct part of aligned_spikes to plot

    xrange = [-.2,+.55];

    for cluster = 1:length(cids)
        % PSTH related variables
        SOM_idx = (stimuli_parameters.Stm.Amplitude == 10); % select stimulus trials
        ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
        preT  = -0.2;
        postT = 0.55;
        binsize = 0.01;

        fig = subplot(2,1,1); % rasterplot
        % set(gcf,'position',[500,150,900,700])
        set(gcf,'position',[500,150,600,400])

        % make rasterplot
        Var = stimuli_parameters.Stm.Amplitude;
        [f, YTick, YTickLab] = plotraster(fig, aligned_spikes(61:140, cluster), Var, [0, 0, 0], [], 1);
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
        % TO DO: histcounts select correct trials
        fig = subplot(2,1,2);

        % index into correct range
        [N, edges] = histcounts(vertcat(aligned_spikes{61:140, cluster}), preT:binsize:postT);
        histogram('BinEdges', edges, 'BinCounts', N, 'FaceColor', '#D95319') % spike/s

        %[N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
        % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
        % hold on
        % [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
        %histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')

        %format axis
        legend('stimulus', 'control')
        xlabel('Time (s)')
        ylabel('Spike rate (Hz)')
        fig.FontSize = 11;
        xlim(fig,xrange);

        sgtitle(['Cluster ' num2str(cids(cluster))])
        %title([num2str(cids(cluster)),' (', num2str(length(spiketimes{cluster})),' total spikes)'])
        % % save plot
        % %figname = sprintf('PSTH_S06_cluster_%i', cluster);
        % %saveas(gcf, fullfile(OutPath, figname));

        disp(stimuli_parameters.Par.Rec)
        disp(NStim)
        disp(stim_counter)
        disp(last_stim)

        % display untill button press
        waitforbuttonpress
        close

    end

%% SOM OLD
    xrange = [-.2,+.55];

    for cluster = 1:length(cids)
        % PSTH related variables
        SOM_idx = (stimuli_parameters.Stm.Amplitude == 10); % select stimulus trials
        ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
        preT  = -0.2;
        postT = 0.55;
        binsize = 0.01;

        fig = subplot(2,1,1); % rasterplot
        % set(gcf,'position',[500,150,900,700])
        set(gcf,'position',[500,150,600,400])

        % make rasterplot
        Var = stimuli_parameters.Stm.Amplitude;
        [f, YTick, YTickLab] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
        yticks(YTick{1});
        yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
        ylim(f,yrange);
        yticklabels(unique(stimuli_parameters.Stm.Amplitude));
        xlim(f,xrange);

        % format axis
        xlabel('Time (s)')
        ylabel('Stimulus off / on')
        fig.FontSize = 11;

        % make histogram / PSTH
        fig = subplot(2,1,2);
        [N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
        histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
        hold on
        [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
        histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')

        %format axis
        legend('stimulus', 'control')
        xlabel('Time (s)')
        ylabel('Spike rate (Hz)')
        fig.FontSize = 11;
        xlim(fig,xrange);

        sgtitle(['Cluster ' num2str(cids(cluster))])
        %title([num2str(cids(cluster)),' (', num2str(length(spiketimes{cluster})),' total spikes)'])
        % % save plot
        % %figname = sprintf('PSTH_S06_cluster_%i', cluster);
        % %saveas(gcf, fullfile(OutPath, figname));

        % display untill button press
        waitforbuttonpress
        close


    end

   