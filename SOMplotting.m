function fig = SOMplotting(sessions, cids, OutPath, BehaviorPath, saveplots)

%load correct files
session = sessions(1); % select experimental session

sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters_som = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes_som = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
aligned_spikes_som = aligned_spikes_som.SpkT;

session = sessions(2); % select control session

sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters_ctrl = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes_ctrl = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
aligned_spikes_ctrl = aligned_spikes_ctrl.SpkT;

%format data to plot together
stimuli_parameters = vertcat(stimuli_parameters_som.Stm, stimuli_parameters_ctrl.Stm);
aligned_spikes = vertcat(aligned_spikes_som, aligned_spikes_ctrl);

xrange = [-.202, +.702];

% plot SOM data based on 'Location'
for cluster = 1:length(cids)

    %figure;
    fig = figure;
    ax = gca;
    %fig = subplot(2,1,1); % rasterplot
    %set(fig,'position',[500,150,800,600]) %set(gcf,'position',[500,150,600,400])
   % location = ~contains(stimuli_parameters.Location, {'air', 'ai1'}); % sort exp session first

    % define stimulus variable space
    if max(stimuli_parameters.SomFreq) == 0
        Var = stimuli_parameters.Amplitude;

        % amplitude options 0 0.05 0.1 0.3 (0.5)
        SOM_on_idx = (stimuli_parameters_som.Stm.Amplitude == 0.3);
        SOM_off_idx = (stimuli_parameters_som.Stm.Amplitude == 0);
        air_on_idx = (stimuli_parameters_ctrl.Stm.Amplitude == 0.3);
        air_off_idx = (stimuli_parameters_ctrl.Stm.Amplitude == 0);
    else
        Var = [stimuli_parameters.SomFreq, stimuli_parameters.Amplitude];
       
        % amplitude options 0 0.1 0.3
        SOM_on_idx = (stimuli_parameters_som.Stm.Amplitude == 0.3);
        SOM_off_idx = (stimuli_parameters_som.Stm.Amplitude == 0);
        air_on_idx = (stimuli_parameters_ctrl.Stm.Amplitude == 0.3);
        air_off_idx = (stimuli_parameters_ctrl.Stm.Amplitude == 0);
    end
    % if strcmp(Par.SomatosensoryActuator, 'Piezo')
    % Var = [stimuli_parameters.Amplitude, location];
    %Var =  [stimuli_parameters.Stm.Freq,stimuli_parameters.Stm.Intensity];

    % [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
    [f, YTick, YTickLab,~,~,YTickLim] = plotraster(ax, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
    yticklabels({'air stim off' 'stim off' 'air stim on' 'stim on'});

    % elseif strcmp(Par.SomatosensoryActuator, 'DC Motor')
    %     Var = [stimuli_parameters.Amplitude, stimuli_parameters.SomFreq, location];
    %     %Var =  [stimuli_parameters.Stm.Freq,stimuli_parameters.Stm.Intensity];
    % 
    %     % [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
    %     [f, YTick, YTickLab,~,~,YTickLim] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [15 30], 1);
    %     yticklabels({'air stim off' 'stim off' 'air stim on' 'stim on'});
    % 
    % end

    yticks(YTick{2});
    yrange = [min(YTick{end}) - 25, max(YTick{end}) + 25];
    ylim(f, yrange);
    xlim(f, xrange);

    % add demarcation lines between categories
    horizontalLine(YTickLim, ax)

    % format axis
    xlabel('Time (s)')
    %ylabel('Actuator position')
    ax.FontSize = 11;
    %title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

    % make histogram / PSTH
    % fig = subplot(2,1,2);
    % preT  = -0.2;
    % postT = 0.7;
    % binsize = 0.01; %0.03; % bit of smoothing psth
    % 
    % % SOM on aligned_spikes_som
    % [N, edges] = histcounts(vertcat(aligned_spikes_som{SOM_on_idx, cluster}), preT:binsize:postT);
    % % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
    % plot(edges(1:end-1),((N/sum(SOM_on_idx == 1))/binsize),'Color', '#00636C','LineWidth',1.5)
    % hold on
    % 
    % % SOM off aligned_spikes_som
    % [N,edges] = histcounts(vertcat(aligned_spikes_som{SOM_off_idx, cluster}), preT:binsize:postT);
    % % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')
    % plot(edges(1:end-1),((N/sum(SOM_on_idx == 0))/binsize),'Color', '#4FC5D3','LineWidth',1.5)
    % 
    % % air on aligned_spikes_ctrl
    % [N, edges] = histcounts(vertcat(aligned_spikes_ctrl{air_on_idx, cluster}), preT:binsize:postT);
    % % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
    % plot(edges(1:end-1),((N/sum(air_on_idx == 1))/binsize),'Color', '#501574','LineWidth',1.5)
    % 
    % % air off
    % [N, edges] = histcounts(vertcat(aligned_spikes_ctrl{air_off_idx, cluster}), preT:binsize:postT);
    % % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
    % plot(edges(1:end-1),((N/sum(air_on_idx == 1))/binsize),'Color', '#B673C8','LineWidth',1.5)
    % 
    % %format axis
    % legend('stimulus on', 'stimulus off', 'air stimulus on', 'air stimulus off')
    % xlabel('Time (s)')
    % ylabel('Spike rate (Hz)')
    % fig.FontSize = 11;
    % xlim(fig,xrange);

    %sgtitle(['                    Fur stroke (unit ' num2str(cids(cluster)) ')'])
    sgtitle([stimuli_parameters.Location(1) 'unit' num2str(cids(cluster))])

    % save plot
    if saveplots
        figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', ...
            str2double(stimuli_parameters_som.Par.MouseNum), str2double(stimuli_parameters_som.Par.Set), stimuli_parameters_som.Par.Rec, cids(cluster));
        saveas(gcf, fullfile(OutPath, [figname '.jpg']));
        saveas(fig, fullfile(OutPath, figname));
    end

    hold off

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
