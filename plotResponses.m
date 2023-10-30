function fig = plotResponses(BehaviorPath, relevant_sessions, aligned_spikes, cids, plotSession)
% plot data: rater & PSTH

% check input
if isempty(aligned_spikes) % to add: look for saved file in folder
    error('no file found')
    %dir(fullfile(uigetdir('D:\')));
    % elseif BehaviorPath ~= 7
    %     error('Folder to stimulus files not found.')
elseif isempty(relevant_sessions)
    relevant_sessions = input('Enter first and last session number in this recording: ');
end

stim_files = dir(fullfile(BehaviorPath, '\*.mat'));
stim_counter = 1;

% select correct stimuli & response group (per session)
for file = relevant_sessions(1):relevant_sessions(2)

    % load correct trial and initiate stimuli counter
    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);
    NStim = size(stimuli_parameters.Stm, 1); % nr trials to align per stim session
    last_stim = stim_counter + NStim - 1;

    % plot FRA data
    if strcmp(stimuli_parameters.Par.Rec, 'FRA') && contains(plotSession, 'FRA')
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
    end

    % plot SOM data
    if strcmp(stimuli_parameters.Par.Rec, 'SOM') && contains(plotSession, 'SOM')

        xrange = [-.2, +.55];

        for cluster = 1:length(cids)
         
            fig = subplot(2,1,1); % rasterplot
            % set(gcf,'position',[500,150,900,700])
            set(gcf,'position',[500,150,600,400])

            % make rasterplot
            Var = stimuli_parameters.Stm.Amplitude;
            %[f, YTick, YTickLab] = plotraster(fig, aligned_spikes(1:40, 1), Var, [0, 0, 0], [], 1);
            [f, YTick, YTickLab] = plotraster(fig, aligned_spikes(stim_counter:last_stim, cluster), Var, [0, 0, 0], [], 1);
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
            % PSTH related variables
            SOM_idx = (stimuli_parameters.Stm.Amplitude == 10); % select stimulus trials
            ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
            preT  = -0.2;
            postT = 0.55;
            binsize = 0.01;

            fig = subplot(2,1,2);

            % index into correct range
            %  [N, edges] = histcounts(vertcat(aligned_spikes{1:40, 1}), preT:binsize:postT);
            [N, edges] = histcounts(vertcat(aligned_spikes{stim_counter:last_stim, cluster}), preT:binsize:postT);
            histogram('BinEdges', edges, 'BinCounts', N, 'FaceColor', '#D95319') % spike/s

            % [N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
            % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
            % hold on
            % [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
            % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')

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
            waitforbuttonpress
            close

        end
    end

    % plot noise data
    if strcmp(stimuli_parameters.Par.Rec, 'AMn') && contains(plotSession, 'noise')

        xrange = [-.2, +.4];

        for cluster = 1:length(cids)

            fig = subplot(2,1,1); % rasterplot
            %set(gcf,'position',[500,150,900,700])
            set(gcf,'position',[500,150,600,400])

            % make rasterplot
            Var = stimuli_parameters.Stm.Intensity;
            [f, YTick, YTickLab] = plotraster(fig, aligned_spikes(stim_counter:last_stim, cluster), Var, [0, 0, 0], [], 1);
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

            % set range & figure
            % preT  = -0.2;
            % postT = 0.4;
            % binsize = 0.01;

            fig = subplot(2,1,2);

            [N, edges] = histcounts(vertcat(aligned_spikes{stim_counter:last_stim, cluster}), preT:binsize:postT);
            histogram('BinEdges', edges, 'BinCounts', N, 'FaceColor', '#D95319') % spike/s

            % display untill button press
            waitforbuttonpress
            close

        end
    end

    stim_counter = stim_counter + NStim; % keep track of how many stimuliu have past

end
end
