function fig = SOMplotting(stimuli_parameters, aligned_spikes, cids, OutPath, saveplots)

% define shared x-lim parameters
preT  = -0.2;
postT = 0.7;
xrange = [preT, postT];
binsize = 0.01;

start_stim = 0;
end_stim = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
xlinerange = [start_stim end_stim];


% plot SOM data based on 'Location'
for cluster = 1:length(cids)

    %fig = figure;
    %ax = gca;
    %set(fig,'position',[500,150,800,600])

    % select 1 frequency to plot
    all_freqs = unique(stimuli_parameters.Stm.SomFreq);
    for freq = 1:length(all_freqs)

        % initiate figure
        fig = subplot(2,1,1); % rasterplot

        % define stimulus variable space
        if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
            Var = stimuli_parameters.Stm.Amplitude;
            yaxislabels = unique(stimuli_parameters.Stm.Amplitude);
            yaxistext = 'Pressure (V)';
        elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
            index = stimuli_parameters.Stm.SomFreq == all_freqs(freq);
            SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
            SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
            Var = [SOM_Hz, SOM_Amp];
            %Var = [SOM_Hz, stimuli_parameters.Stm.Amplitude];
            yaxislabels = unique(stimuli_parameters.Stm.Amplitude);
            %yaxistext = 'Vibrotactile stimulation (Hz)';
            yaxistext = [num2str(all_freqs(freq)) ' Hz'];
        end

        % [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(fig, aligned_spikes(index, cluster), Var, [0, 0, 0], [], 1);
        yticks(YTick{2});
        yrange = [min(YTick{end}) - 15, max(YTick{end}) + 15];
        ylim(f,yrange);
        xlim(f,xrange);

        % add demarcation lines
        xline(xlinerange) % on/off set
        horizontalLine(YTickLim,fig) % between categories

        % format axis
        xlabel('Time (s)')
        yticklabels(yaxislabels)
        ylabel(yaxistext)
        %fig.FontSize = 11;
        title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        % make histogram / PSTH
        fig = subplot(2,1,2);

        % plot low and high amplitude seperate
        SOM_Amp_low = (stimuli_parameters.Stm.Amplitude == 0.1) & (stimuli_parameters.Stm.SomFreq == all_freqs(freq));
        SOM_Amp_high = (stimuli_parameters.Stm.Amplitude == 0.3) & (stimuli_parameters.Stm.SomFreq == all_freqs(freq));

        % [N, edges] = histcounts(vertcat(aligned_spikes_som{SOM_on_idx, cluster}), preT:binsize:postT);
        [N, edges] = histcounts(vertcat(aligned_spikes{SOM_Amp_low, cluster}), preT:binsize:postT);
        % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
        % plot(edges(1:end-1),((N/sum(SOM_on_idx == 1))/binsize),'Color', '#00636C','LineWidth',1.5)
        plot(edges(1:end-1),((N/sum(SOM_Amp_low))/binsize),'Color', '#00636C','LineWidth',1)

        hold on
        [N, edges] = histcounts(vertcat(aligned_spikes{SOM_Amp_high, cluster}), preT:binsize:postT);
        plot(edges(1:end-1),((N/sum(SOM_Amp_high))/binsize),'Color', '#D95319','LineWidth',1)
        hold off

        %format axis
        %legend('0.1V', '0.3V', 'Location','northwest')
        xlabel('Time (s)')
        ylabel('Spike rate (Hz)')
        fig.FontSize = 11;
        xlim(fig,xrange);
        %ylim(fig, [0, 300])
        xline(xlinerange) % on/off set


        %sgtitle(['                    Fur stroke (unit ' num2str(cids(cluster)) ')'])
        %sgtitle([stimuli_parameters.Par.SomatosensoryLocation(1) 'unit' num2str(cids(cluster))])

        % save plot
        if saveplots
            figname = sprintf('M%.2i_S%.2i_%s_cluster_%i_%iHz', ...
                str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster), all_freqs(freq));
            saveas(gcf, fullfile(OutPath, [figname '.jpg']));
            %saveas(fig, fullfile(OutPath, figname));
        end

    end

end


function horizontalLine(YTickLim, fig)

for i = 1:size(YTickLim,1)
    yline(fig,YTickLim(i,1)-3,':k');
    yline(fig,YTickLim(i,2)+3,':k');
end

% box
x = -0.2;w=1;
for i = 1:size(YTickLim,1)
    r(i) = rectangle(fig,'Position',[x,YTickLim(i,1)-0.5,w,YTickLim(i,2)-YTickLim(i,1)+1],...
        'FaceColor',[0,0,0,0],'EdgeColor','none');
end
