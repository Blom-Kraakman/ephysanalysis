function fig = SOMplotting(stimuli_parameters, aligned_spikes, cids, OutPath, saveplots, condition)

% define shared x-lim parameters
preT  = -0.2;
postT = 0.7;
xrange = [preT, postT];
binsize = 0.01; %10ms

start_stim = 0;
end_stim = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
xlinerange = [start_stim end_stim]; %end_stim
%start_som = max(stimuli_parameters.Stm.SomAudSOA)/1000;


% plot SOM data based on 'Location'
for cluster = 1:length(cids)

    % % add delay to "SO","OO" if needed
    % idx = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
    % for ii = idx'
    %     aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} + start_som;
    % end

    if cids(cluster)

        %fig = figure;
        %ax = gca;
        %set(fig,'position',[500,150,800,600])

    % select 1 frequency to plot
    all_freqs = unique(stimuli_parameters.Stm.SomFreq);
    for freq = 1:length(all_freqs)

        % initiate figure
        %fig = figure('position',[-918,294,520,370]);
        fig = figure;
        raster_fig = subplot(4,1,4); % rasterplot

        % define stimulus variable space
        if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
            index = ones(size(stimuli_parameters.Stm,1),1);
            Var = stimuli_parameters.Stm.Amplitude;
            yaxislabels = unique(stimuli_parameters.Stm.Amplitude);
            yaxistext = 'Pressure (V)';
            linecolor = 'k';
        elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
            % define condition and shift alignspikes if needed
            if strcmp(condition, 'SA')
                index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.3) & (stimuli_parameters.Stm.AudIntensity ==45); % select only high pressure case
                linecolor = '#742f9e';

                SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
                SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
                Var = [SOM_Hz, SOM_Amp];
            elseif strcmp(condition, 'SO')
               % index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.3) & (stimuli_parameters.Stm.AudIntensity == -Inf); % select only high pressure case

                index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.3170) & (stimuli_parameters.Stm.AudIntensity == -Inf); % select only high pressure case
                linecolor = '#52A3CF';

                SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
                SOM_Amp = stimuli_parameters.Stm.Amplitude(index);
                Var = [SOM_Hz, SOM_Amp];
            elseif strcmp(condition, 'OA')
                index = strcmp(stimuli_parameters.Stm.MMType, 'OA'); % select only high pressure case
                linecolor = '#B61768';

                %SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
                %SOM_Amp = stimuli_parameters.Stm.Amplitude(index);

                Var = [stimuli_parameters.Stm.AudIntensity(index)];
            end

            
            %Var = [SOM_Hz, stimuli_parameters.Stm.Amplitude];
            yaxislabels = unique(stimuli_parameters.Stm.Amplitude);
            %yaxistext = 'Vibrotactile stimulation (Hz)';
            yaxistext = [num2str(all_freqs(freq)) ' Hz'];

        end

        if ~sum(index)
            close(gcf)
            continue
        end

        % [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(raster_fig, aligned_spikes(index, cluster), Var, [0, 0, 0], [], 1);
        %yticks(YTick{2});
        yticks([1 20]);
        yrange = [min(YTick{end}) - 1, max(YTick{end}) + 1]; % 5 ipv 15
        ylim(f,yrange);
        xlim(f,xrange);

        % add demarcation lines
        xline(xlinerange) % on/off set
        horizontalLine(YTickLim,raster_fig) % between categories

        % format axis
        xlabel('Time (s)')
        ylabel('trials')
        raster_fig.FontSize = 11;
        %sgtitle(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

        % make histogram / PSTH
        psth_fig = subplot(4,1,2:3);

        % plot low and high amplitude seperate
        %SOM_Amp_low = (stimuli_parameters.Stm.Amplitude == 0.1) & (stimuli_parameters.Stm.SomFreq == all_freqs(freq));
        SOM_Amp_high = index; %(stimuli_parameters.Stm.Amplitude == 0.3) & (stimuli_parameters.Stm.SomFreq == all_freqs(freq));

        % [N, edges] = histcounts(vertcat(aligned_spikes_som{SOM_on_idx, cluster}), preT:binsize:postT);
        % plot(edges(1:end-1),((N/sum(SOM_on_idx == 1))/binsize),'Color', '#00636C','LineWidth',1.5)

        %[N, edges] = histcounts(vertcat(aligned_spikes{SOM_Amp_low, cluster}), preT:binsize:postT);
        % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
        %plot(edges(1:end-1),((N/sum(SOM_Amp_low))/binsize),'Color', '#00636C','LineWidth',1)

        hold on
        [N, edges] = histcounts(vertcat(aligned_spikes{SOM_Amp_high, cluster}), preT:binsize:postT);
        plot(edges(1:end-1),((N/sum(SOM_Amp_high))/binsize),'Color', linecolor,'LineWidth',2)
        hold off

        %format axis
        %legend('0.1V', '0.3V', 'Location','northwest')
        %xlabel('Time (s)')
        ylabel('Spike rate (Hz)')
        psth_fig.FontSize = 11;
        xlim(psth_fig,xrange);
        ylim(psth_fig, [0, 320])
        xline(xlinerange) % on/off set
        xticklabels(psth_fig, {[]})

        % plot sine wave
        subplot(4,1,1)
        Amplitude = 5;
        %StimOnset = 250;
        [som_waveform,tt_stim] = gensomwaveform('UniSine',500, Amplitude,all_freqs(freq),0,30000);
        plot(tt_stim+0.250,som_waveform, 'k', 'LineWidth',1); %hold(ax(1),'on');
        % plot(tt_stim+0.250,som_waveform, 'k', 'LineWidth',1); %hold(ax(1),'on');

        %xlim()
        xlim(xrange)
        axis off

        sgtitle([num2str(all_freqs(freq)) ' Hz']);

        % save plot
        if saveplots
            figname = sprintf('M%.2i_S%.2i_%s_cluster_%i_%iHz_%s', ...
                str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster), all_freqs(freq), condition);
            saveas(gcf, fullfile(OutPath, [figname '.jpg']));
            saveas(fig, fullfile(OutPath, figname));
        end
    end

    end

%close all
end
end

%% local functions

% copied from NeuroPassiveFcns.m
function [som_waveform,tt] = gensomwaveform(Waveform,StimDur, Amplitude,SomFreq,Ramp,Fs)
    StimDurSamp = ceil(StimDur * 0.001 * Fs);
    som_waveform = nan(1,StimDurSamp);
    tt = 0:1/Fs:(StimDur* 0.001);
    switch Waveform
        case {'Square'}
            som_waveform = Amplitude .* ones(1,StimDurSamp);
            som_waveform(1) = 0; som_waveform(end) = 0; % zero at beginning or end
        case {'UniSine'}
            som_waveform =  0.5 * Amplitude .* ( 1-cos(2*pi*SomFreq .* tt) );
        case {'BiSine'}
            som_waveform =  Amplitude .* ( sin(2*pi*SomFreq .* tt) );
    end

    % apply envelope (On-/Off-ramps)
    if Ramp > 0
        Nenv			=	round( Ramp * 10^-3 * Fs );
	    som_waveform    =	envelope(som_waveform',Nenv)';
    end

end

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
