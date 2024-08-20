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

% plot noise and AM
if strcmp(stimuli_parameters.Par.Rec, 'AMn')

    % define shared x-lim parameters
    preT  = -0.2;
    postT = 1.4;
    xrange = [preT, postT];
    binsize = 0.01;
    xlinerange = [0 0.5 1];
    % str2double(stimuli_parameters.Par.AMStimTime)/1000

    for cluster = 1:length(cids)

        fig = figure;
        ax = gca;

        % fig = subplot(2,1,1); % rasterplot
        set(gcf,'position',[500,150,900,700])

        if max(stimuli_parameters.Stm.Mf) > 0
            Var = [stimuli_parameters.Stm.Intensity, stimuli_parameters.Stm.Mf];
            ytick_label = unique(stimuli_parameters.Stm.Mf);
            y_tick_level = 2;
        elseif max(stimuli_parameters.Stm.Mf) == 0
            Var = stimuli_parameters.Stm.Intensity;
            ytick_label = unique(stimuli_parameters.Stm.Intensity);
            y_tick_level = 1;
        end

        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(ax, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
        yticks(YTick{y_tick_level});
        yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
        ylim(f,yrange);
        yticklabels(ytick_label);
        xlim(f,xrange);

        % add demarcation lines
        xline(xlinerange) % on/off set
        horizontalLine(YTickLim, ax) % between categories

        % format axis
        xlabel('Time (s)')
        ylabel('AM frequency (Hz)')
        %ax.FontSize = 11;
        %title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        % if max(stimuli_parameters.Stm.Mf) > 0
        %     %     AM_0_idx = stimuli_parameters.Stm.Mf == 0;
        %     %     AM_idx = stimuli_parameters.Stm.Mf ~= 0;
        %     %
        %     %     % set range & figure
        %     %     fig = subplot(2,1,2);
        %     %
        %     %     % make histogram / PSTH
        %     %     [N, edges] = histcounts(vertcat(aligned_spikes{AM_0_idx, cluster}), preT:binsize:postT);
        %     %     plot(edges(1:end-1), ((N/sum(AM_0_idx))/binsize), 'Color', '#0072BD', 'LineWidth',1.5) % spike/s
        %     %     hold on
        %     %
        %     %     [N,edges] = histcounts(vertcat(aligned_spikes{(AM_idx), cluster}), preT:binsize:postT);
        %     %     plot(edges(1:end-1), ((N/sum(AM_idx))/binsize), 'Color', '#D95319', 'LineWidth',1.5) % spike/s
        %     %
        %     %     %format axis
        %     %     legend('No modulation', 'Amplitude modulation')
        % 
        % else
        %     % set variables for PSTH
        %     AM_15_idx = stimuli_parameters.Stm.Intensity == 15; % select trials
        %     AM_30_idx = stimuli_parameters.Stm.Intensity == 30;
        %     AM_45_idx = stimuli_parameters.Stm.Intensity == 45;
        %     AM_60_idx = stimuli_parameters.Stm.Intensity == 60;
        %     AM_0_idx = (stimuli_parameters.Stm.Intensity ~= 15) & (stimuli_parameters.Stm.Intensity ~= 30) & ...
        %         (stimuli_parameters.Stm.Intensity ~= 45) & (stimuli_parameters.Stm.Intensity ~= 60);
        % 
        %     % set range & figure
        %     fig = subplot(2,1,2);
        % 
        %     % make histogram / PSTH
        %     [N, edges] = histcounts(vertcat(aligned_spikes{AM_60_idx, cluster}), preT:binsize:postT);
        %     plot(edges(1:end-1), ((N/sum(AM_60_idx))/binsize), 'Color', '#0072BD', 'LineWidth',1.5) % spike/s
        %     %plot(edges(1:end-1), ((N/sum(AM_60_idx))/binsize), 'Color', '#4FC5D3', 'LineWidth',1.5) % spike/s
        %     hold on
        % 
        %     [N,edges] = histcounts(vertcat(aligned_spikes{(AM_45_idx), cluster}), preT:binsize:postT);
        %     plot(edges(1:end-1), ((N/sum(AM_45_idx))/binsize), 'Color', '#D95319', 'LineWidth',1.5) % spike/s
        %     %plot(edges(1:end-1), ((N/sum(AM_45_idx))/binsize), 'Color', '#B673C8', 'LineWidth',1.5) % spike/s
        % 
        %     [N,edges] = histcounts(vertcat(aligned_spikes{(AM_30_idx), cluster}), preT:binsize:postT);
        %     plot(edges(1:end-1), ((N/sum(AM_30_idx))/binsize), 'Color', '#7E2F8E', 'LineWidth',1.5) % spike/s
        % 
        %     [N,edges] = histcounts(vertcat(aligned_spikes{(AM_15_idx), cluster}), preT:binsize:postT);
        %     plot(edges(1:end-1), ((N/sum(AM_15_idx))/binsize), 'Color', '#77AC30', 'LineWidth',1.5) % spike/s
        % 
        %     [N, edges] = histcounts(vertcat(aligned_spikes{AM_0_idx, cluster}), preT:binsize:postT);
        %     plot(edges(1:end-1), ((N/sum(AM_0_idx))/binsize), 'Color', 'k', 'LineWidth',1.5) % spike/s
        % 
        %     %format axis
        %     legend('60dB (SPL)', '45dB (SPL)', '30dB (SPL)', '15dB (SPL)', 'no sound')
        % 
        % end

        %xlabel('Time (s)')
        %ylabel('Spike rate (Hz)')
        %xlim(fig,xrange);
        sgtitle(['Broadband noise (unit ' num2str(cids(cluster)) ')'])
        figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        saveas(fig, fullfile(OutPath, figname));

        hold off

    end
end

% plot SOM data
if strcmp(stimuli_parameters.Par.Rec, 'SOM')

    % define shared x-lim parameters
    preT  = -0.2;
    postT = 0.7;
    xrange = [preT, postT];
    binsize = 0.01;

    start_stim = 0;
    end_stim = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
    xlinerange = [start_stim end_stim];
    
    for cluster = 1:length(cids)

        fig = figure;
        ax = gca;
        %fig = subplot(2,1,1); % rasterplot
        %set(gcf,'position',[500,150,800,600])

        % define stimulus variable space
        if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
            Var = stimuli_parameters.Stm.Amplitude;
            yaxislabels = unique(stimuli_parameters.Stm.Amplitude);
            yaxistext = 'Pressure (V)';
        elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
            Var = [stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude];
            yaxislabels = unique(stimuli_parameters.Stm.SomFreq);
            yaxistext = 'Vibrotactile stimulation (Hz)';
        end

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
if strcmp(stimuli_parameters.Par.Rec, 'SxA')

    % define shared x-lim parameters
    preT  =  max(-str2double(stimuli_parameters.Par.SomatosensoryISI)/2000, -0.2);
    postT = min((max(str2double(stimuli_parameters.Par.AuditoryStimTime), str2double(stimuli_parameters.Par.SomatosensoryStimTime)) ...
            + str2double(stimuli_parameters.Par.SomatosensoryISI)/2)/1000, 1.5);
    xrange = [preT, postT];
    binsize = 0.01;
    start_aud = 0;
    start_som = max(stimuli_parameters.Stm.SomAudSOA)/1000;
    end_som = max(stimuli_parameters.Stm.SomAudSOA)/1000 + str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
    end_aud = str2double(stimuli_parameters.Par.AuditoryStimTime)/1000;

    % start_aud = 0;
    % start_som = 0.25;
    % end_som = 0.75;
    % end_aud = 1;

    xlinerange = [start_aud start_som end_som end_aud];

    for cluster = 1:length(cids)

        fig = figure;
        ax = gca;

        %fig = subplot(2,1,1); % rasterplot
        set(gcf,'position',[500,150,900,700])

        % make rasterplot
        raster_color = [0, 0, 0];
        raster_yinc = [];%1,1,1];
        % define stimulus variable space
        if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
            Var = [stimuli_parameters.Stm.Var25, stimuli_parameters.Stm.AudIntensity, stimuli_parameters.Stm.Amplitude];
            yaxislabels = {'Control', 'Pressure only', 'Pressure + noise' 'Noise only'};
            yaxistext = 'Pressure (V)';
            raster_yinc = [10,20,20];

            % add delay to "SO","OO" if needed
            idx = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
            for ii = idx'
                aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} + start_som;
            end
        elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') || strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'BiSine')
            % Var = [stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude ~=0.1,stimuli_parameters.Stm.Var25];
            Var = [stimuli_parameters.Stm.Var25, stimuli_parameters.Stm.SomFreq, stimuli_parameters.Stm.Amplitude];
            raster_yinc = [5,10,10];
            % yaxislabels = unique(stimuli_parameters.Stm.SomFreq);
            yaxislabels = {'Control', 'Vibrotactile only', 'Vibrotactile + noise' 'Noise only'};
            %yaxistext = 'Vibrotactile stimulation (Hz)';

            % add delay to "SO","OO" if needed
            idx = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
            for ii = idx'
                aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} + start_som;
            end
        end

        % [f, YTick, YTickLab] = plotraster(gca, aligned_spikes(:, 1), Var, [0, 0, 0], [], 1);
        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(ax, aligned_spikes(:, cluster), Var, raster_color, raster_yinc, 1);
        yticks(YTick{1});
        yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];        ylim(f,yrange);
        xlim(f,xrange);

        % add demarcation lines
        xline(xlinerange) % on/off set
        horizontalLine(YTickLim, ax) % between categories

        % format axis
        xlabel('Time (s)')
        yticklabels(yaxislabels)
        %ylabel(yaxistext)
        %fig.FontSize = 11;
        title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        %save figure
        figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        saveas(fig, fullfile(OutPath, figname));
    end

end
end

function horizontalLine(YTickLim, fig)

for i = 1:size(YTickLim,1)
    yline(fig,YTickLim(i,1)-3,':k');
    yline(fig,YTickLim(i,2)+3,':k');
end

% % box
% x = -0.2; w=1;
% for i = 1:size(YTickLim,1)
%     %r(i) = rectangle(fig,'Position',[x,YTickLim(i,1)-0.5,w,YTickLim(i,2)-YTickLim(i,1)+1],...
%      %   'FaceColor',[0,0,0,0],'EdgeColor','none');
% end
end