% channel map
%xcoords = channel_positions(:,1);
%ycoords = channel_positions(:,2);
figure
margin = 50;
xMin = min(xcoords); xMax = max(xcoords);
xSpan = xMax-xMin; xMid = 0.5*(xMax+xMin);
yMin = min(ycoords); yMax = max(ycoords);
ySpan = yMax-yMin; yMid = 0.5*(yMax+yMin);
maxSpan = max(xSpan,ySpan);

figa = scatter(xcoords, ycoords, ".", 'k');
hold on;
text(xcoords, ycoords,num2str(chanMap))

%% Parameters
f = 10;        % Frequency in Hz
Fs = 1000;     % Sampling frequency in Hz
T = 0.5;       % Duration in seconds
t = -0.2:1/Fs:0.7;  % Time vector from -0.2 to 0.7 seconds

% Sine wave
y = sin(2 * pi * f * t);

% Unidirectional sine wave (offset by 1)
y_unidirectional = (y + abs(min(y))) .* (t >= 0 & t <= T);

% Adding stable line at 0 between -0.2 and 0 seconds, and between 0.5 and 0.7 seconds
y_unidirectional(t < 0 | t > T) = 0;

% Plotting
figure;
plot(t, y_unidirectional);


%% define start/end cycle + cycle window
% does not work well

nCycles = (stimuli_parameters.Stm.SomFreq .* (str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000));
cycle_window = ((str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000) ./ nCycles); %samples
trialDur = ((str2double(stimuli_parameters.Par.SomatosensoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI))/1000) * Fs; %in samples

%cycle_aligned = nan(nClusters, nTrials, max(nCycles)); 
cycle_aligned = [];
% get samples corresponding to each cycle window
% for each cell
for cluster = 1:nClusters
    % for all trials
    for stim = 1:nTrials
        % time window each cycle
        cycle_on = 0;
        cycle_off = cycle_window(stim);
        for cycle = 1:nCycles(stim)
            % get spike times in this time window from aligned_spikes
            taligned_spikes = aligned_spikes.SpkT{stim, cluster};
            idx = (taligned_spikes >= cycle_on) & (taligned_spikes <= cycle_off);

            if max(idx) == 0; continue; end

            tcycle_aligned{(cycle)} = taligned_spikes(idx); % store spike times in cell array

            %cycle_aligned(cluster, stim, cycle) = taligned_spikes(idx);

            % update window
            cycle_on = cycle_on + cycle_window(stim);
            cycle_off = cycle_off + cycle_window(stim);
        end

        cycle_aligned = [cycle_aligned, tcycle_aligned];
    end

end


% match matrix with firing rate

%% OLD - remove double spikes from originating in Phy - No longer needed OLD
% make into function
% get all spiketimes from each good unit
spiketimes2 = cell(length(cids), 1);

for cluster = 1:length(cids)
    Tspkt = spiketimes{cluster};
    ISI = zeros(length(Tspkt),1);

    %find minimal time between two spikes
    for spike = 1:(length(Tspkt)-1)
        ISI(spike) = min(abs(Tspkt(spike) - Tspkt(spike+1)));
    end

    % select spikes following next one with more than 0.2ms (6 samples)
    index = ISI>=6;
    spiketimes2{cluster} = Tspkt(index);
    disp(size(spiketimes{cluster}, 1) - size(spiketimes2{cluster}, 1))
end

%spiketimes = spiketimes2;

%% plot aligned spikes raster of each condition
all_freqs = unique(stimuli_parameters.Stm.SomFreq);

for cluster = 1:nClusters
    %for freq = 1:length(all_freqs)
freq = 7;
        figure;

        index = (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == 0.1);
        SOM_Hz = stimuli_parameters.Stm.SomFreq(index);
        Var = SOM_Hz;
        yaxistext = [num2str(all_freqs(freq)) ' Hz'];

        % make rasterplot
        [f, YTick, ~, ~, ~, YTickLim] = plotraster(gca, aligned_spikes(index, cluster), Var, [0, 0, 0], [5 5], 1);
        yrange = [min(YTick{end}) - 50, max(YTick{end}) + 50];

    %end
end

%% use full words TTLs for spike alignment
full_words = readNPY('D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data-A\TTL\full_words.npy');

% gives 1 when either and both stim (2: aud, 16(2^4): som) are on/high
stim_on = bitand(full_words,(2+16)) > 0; 
stim_rise = [true;diff(stim_on) > 0];
stim_fall = [false;diff(stim_on) < 0];
%% check nr ttls in som+aud
TTL_samples = readNPY([TTLPath 'sample_numbers.npy']); % sample nr all recorded TTLs
TTL_states = readNPY([TTLPath 'states.npy']);
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));
stimuli_parameters = load([stim_files(1).folder '\' stim_files(1).name]);

% remove camera TTLs
index = (abs(TTL_states) == 8);
TTL_states(index) = [];
TTL_samples(index) = [];

% remove any non aud or som ttls
index = (abs(TTL_states) ~= 2 & abs(TTL_states) ~= 5);
TTL_states(index) = [];
TTL_samples(index) = [];

%% regular expressions
expression = 'end s';
startpoint = strfind(msgs(1), expression);

extractAfter(msgs(16), expression)


%%test new function
OutPath = 'D:\DATA\Processed\M7';
stim_files = dir(fullfile(OutPath, '*_cluster.fig'));
for filenr = 1:length(stim_files)
    file = openfig([stim_files(filenr).folder '\' stim_files(filenr).name]);
    waitforbuttonpress;
    close(file)
end

%% open multiple figures
[figures, pathname,]=uigetfile('directory','*.fig','Multiselect','on');
for x = 1:length(figures)
    Multi_Figs = [pathname,filesep,figures{x}];
    Op = openfig(Multi_Figs);
end

%% previous version of reading messages from OE
% function sessions_TTLs = makeSessionTTLs(messagesPath)
% import .cvs with text messages
message_samples = readNPY([messagesPath 'sample_numbers.npy']); % session TTLs
Pathmessage_center_text = 'D:\DATA\Processed\M04_message_text.csv';
message_center_text= readtable(Pathmessage_center_text,'ReadVariableNames',false,'Format','%s','Delimiter',',');
%test = table2cell(message_center_text);
sessions_TTLs = table2cell(message_center_text);
for i = 1:length(sessions_TTLs)
    if contains(test{i,1}, 'start')
        sessions_TTLs(i,2) = 1;
    elseif contains(test{i,1}, 'end')
        sessions_TTLs(i,2) = 0;
    else
        error('Error in makeSessionTTLs, check .csv file')
    end
end
% sessions_TTLs(:,1) = [1 1 2 2 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 19 19 20 20 21 21 22 22 23 23]; % session nr
% sessions_TTLs(:,2) = [1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0]; % session start/end (1/0)
sessions_TTLs(:,3) = message_samples(1:length(sessions_TTLs)); % sample nr

%% read .dat files 
%datafile = readtable();
fid         = fopen('D:\DATA\EphysRecordings\M7\M07_2024-02-01_14-09-33\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B.dat', 'r');
fread('D:\DATA\EphysRecordings\M7\M07_2024-02-01_14-09-33\Record Node 108\experiment1\recording1\continuous\Intan-100.Rhythm Data-B.dat', '*int16');
%%
logical(connected(65:72)) = 0;
%%
%chanMapTest = chanMap;

channelstoadd = [65 66 67 68 69 70 71 72]';
kcoordstoadd = NaN(8,1);
xcoordstoadd = [-22.5000 -22.5000 -22.5000 -22.5000 -22.5000 -22.5000 -22.5000 -22.5000]';
ycoordstoadd = [-10 -15 -20 -25 -30 -35 -40 -45]';

kcoords = [kcoords; kcoordstoadd];
xcoords = [xcoords; xcoordstoadd];
ycoords = [ycoords; ycoordstoadd];
chanMap0ind = [chanMap0ind; (channelstoadd-1)];
chanMap = [chanMap; channelstoadd];
connected = logical([connected; 1; 1; 1; 1; 1; 1; 1; 1;]);

save('C:\Users\TDT\Documents\GitHub\KilosortSettings\Kilosort3Config\Intan_cambridgeneurotech_ASSY-77-H5_72ch.mat', "chanMap", "chanMap0ind", "connected", "kcoords", "xcoords", "ycoords", '-mat');
%%
channel_map = readNPY('D:\DATA\EphysRecordingsSorted\M07\channel_map.npy');
channel_positions = readNPY('D:\DATA\EphysRecordingsSorted\M07\channel_positions.npy');

xcoords = channel_positions(:,1);
ycoords = channel_positions(:,2);
margin = 50;
xMin = min(xcoords); xMax = max(xcoords);
xSpan = xMax-xMin; xMid = 0.5*(xMax+xMin);
yMin = min(ycoords); yMax = max(ycoords);
ySpan = yMax-yMin; yMid = 0.5*(yMax+yMin);
maxSpan = max(xSpan,ySpan);

fig = scatter(xcoords, ycoords, ".", 'k');
hold on;

for i = 1:length(channel_map)
    text((channel_positions(i,1)+1), (channel_positions(i,2)+1),num2str(channel_map(i)))
end

% format figure
title('Channel map')
ylabel('Relative depth')
fig.SizeData = 100;
fontsize(13,"points")
axis square
xlim([xMid - 0.5*maxSpan - margin, xMid + 0.5*maxSpan + margin]);
ylim([yMid - 0.5*maxSpan - margin, yMid + 0.5*maxSpan + margin]);
xticks(unique(xcoords));
yticks(unique(ycoords));

hold off;

%%

message_samples = readNPY([messagesPath 'sample_numbers.npy']); % session TTLs
Pathmessage_center_text = 'D:\DATA\Processed\M04_message_text.csv';
message_center_text = readtable(Pathmessage_center_text,'ReadVariableNames',false,'Delimiter',',');
%test = table2cell(message_center_text);
message_center_text = table2cell(message_center_text);
sessions_TTLs = (1:height(message_center_text))';
for i = 1:length(sessions_TTLs)
    % indicate session start (1) and end (0)
    if contains(message_center_text{i,1}, 'start')
        sessions_TTLs(i,2) = 1;
    elseif contains(message_center_text{i,1}, 'end')
        sessions_TTLs(i,2) = 0;
    else
        error('Error in makeSessionTTLs, check .csv file')
    end

    % session number
    for j = 1:27
        if contains(message_center_text{5,1}, 'j')
            sessions_TTLs(:,3) = j;
        end
    end

end

% sessions_TTLs(:,1) = [1 1 2 2 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 19 19 20 20 21 21 22 22 23 23]; % session nr
% sessions_TTLs(:,2) = [1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0]; % session start/end (1/0)
sessions_TTLs(:,4) = message_samples(1:length(sessions_TTLs)); % add sample nr

% manual making of sessions_TTLs
message_samples = readNPY([messagesPath 'sample_numbers.npy']); % session TTLs
sessions_TTLs = [];
sessions_TTLs(:,1) = [1 2 3 3 4 4 5 5 6 6 7 7 8 8 10 10 11 11 12 12 13 13 14 14 17 17]; % session nr
sessions_TTLs(:,2) = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]; % session start/end (1/0)
sessions_TTLs(:,3) = message_samples(1:length(sessions_TTLs)); % sample nr

%% TTLs recorded
% to read text: in python save as csv. perhaps also combine all measures in
% 1 csv for ease later
%regex (regular expression)

% load TTLs + sample numbers
rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
message_samples = readNPY('D:/DATA/EphysRecordings/M4/M04_2023-12-14_12-47-41/Record Node 103/experiment1/recording1/events/MessageCenter/sample_numbers.npy'); % session TTLs
TTL_samples = readNPY([TTLPath 'sample_numbers.npy']); % sample nr all recorded TTLs
TTL_states = readNPY([TTLPath 'states.npy']);
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

%make sessions type + nr trials table
Nr_sessions = relevant_sessions';

for file = 1:length(Nr_sessions)
    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);

    Nr_sessions(file,1) = str2double(stimuli_parameters.Par.Set);
    Nr_sessions(file,2) = size(stimuli_parameters.Stm, 1);

    % notate TTL nr (TTLS: 5 = SOM, 2 = AUD, 8 = CAM)
    if strcmp(stimuli_parameters.Par.Rec, 'FRA')
        Nr_sessions(file,3) = 2;
    elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
        Nr_sessions(file,3) = 2;
    elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
        Nr_sessions(file,3) = 5;
    end

end

% recording sessions table
sessions_TTLs(:,1) = [1 1 2 2 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 19 19 20 20 21 21 22 22 23 23 24 24 25 25 26 26 27 27]; % session nr
sessions_TTLs(:,2) = [1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]; % session start/end (1/0)
sessions_TTLs(:,3) = message_samples; % sample nr

% remove camera TTLs
index = (abs(TTL_states) == 8);
TTL_states(index) = [];
TTL_samples(index) = [];

%keep only TTLs recorded during specific session
for i = 1:length(sessions_TTLs)
    if sessions_TTLs(i,2) == 1
        session_start = sessions_TTLs(i,3); %start
        session_end = sessions_TTLs(i+1,3); %end
        idx = (TTL_samples >= session_start) & (TTL_samples < session_end);
        tTTL_states = TTL_states(idx);
        tTTL_samples = TTL_samples(idx);
    end

    % keep only session specific TTLs
    TTLnr = Nr_sessions(sessions_TTLs(i,1), 3); % get TTLnr corresponding to session
    index = (abs(tTTL_states) == TTLnr);
    ttTTL_states = tTTL_states(index);
    ttTTL_samples = tTTL_samples(index);

    TTLs_sample = [TTLs_sample, tTTL_samples];

    % session always starts with high/positive TTL number
    if ttTTL_states(1) < 0
        ttTTL_samples(1) = [];
        ttTTL_states(1) = [];
    end

    % from session keep only the relevant TTLs (so aud ttl for aud session)
    Srise = ttTTL_samples(ttTTL_states == TTLnr);
    Sfall = ttTTL_samples(ttTTL_states == -TTLnr);

    % remove artefacts where Srise == Sfall
    minDur = 30 ; % samples (= 1ms)
    idx = find ((Sfall - Srise) < minDur);
    Srise(idx) = [];
    Sfall(idx) = [];

end

%% plot data: raster & PSTH
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

%% check for TTL flickers & remove them
if length(unique(tTTL_samples)) ~= length(tTTL_samples)
    disp('removing sample point duplicates in TTLs')

    [uniqueValues, uniqueIndices, indices] = unique(tTTL_samples);

    % Count occurrences of each sample value
    counts = histcounts(indices, 1:numel(uniqueValues)+1);

    % Identify non-duplicate indices
    nonDuplicateIndices = find(counts == 1);

    % Extract only non-duplicate sample values and corresponding TTL states
    ttTTL_samples = uniqueValues(nonDuplicateIndices);
    ttTTL_states = tTTL_states(ismember(indices, nonDuplicateIndices));

end

%TTL_samples(idx)
sessions_TTLs(2,1)

% select on and offset stimuli (5=SOM, 2=AUD)
Srise = TTL_samples((TTL_states == 2) | (TTL_states == 5)); % Srise = column vector of length NStim where the spikes should be aligned
Sfall = TTL_samples((TTL_states == -2) | (TTL_states == -5)); % Sfall = column vector of length nStm indicating end of stimulus

%% plot summary fig - vibrotactile sessions
% boxplot, per unit

close all
amp = 0.3; 

uFreq = unique(stimuli_parameters.Stm.SomFreq);
nFreq = length(uFreq);
nClusters = length(cids);
condition = nan(nClusters, nFreq);
figure

for cluster = 1:nClusters
    for freq = 1:nFreq
        condition(cluster, freq) = median(stimulusRate([stimuli_parameters.Stm.SomFreq] == uFreq(freq) & [stimuli_parameters.Stm.Amplitude] == amp, cluster))';
    end

    plot(uFreq, mean(condition(cluster, :) - condition(cluster, 1), 1), ':o')
    hold on

end

xticks(uFreq)
xticklabels(uFreq(1:nFreq));
xlabel(['Vibrotactile stimulus (Hz), ' num2str(amp) ' (V)'])
ylabel('Normalized firing rate (Hz)')

title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

save figure
figname = sprintf('M%.2i_S%.2i_%s_mean', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

% plot summary fig: boxplot, all units
close all
condition = nan(nClusters, nFreq);
figure

% select data
for freq = 1:nFreq
    condition(:, freq) = median(stimulusRate([stimuli_parameters.Stm.SomFreq] == uFreq(freq) & [stimuli_parameters.Stm.Amplitude] == amp, :))';
end

% plot
figure
boxplot(condition)

xticklabels(uFreq(1:nFreq));
xlabel(['Vibrotactile stimulus (Hz), ' num2str(amp) ' (V)'])
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

% save figure
figname = sprintf('M%.2i_S%.2i_%s_boxplot_all units', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));


%% plot summary fig - pressure sessions
% boxplot, per  units
close all

% parameters
uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
nClusters = length(cids);
condition = nan(nClusters, nAmp);

figure
for cluster = 1:nClusters

    % select data
    for amplitude = 1:nAmp
        condition(cluster, amplitude) = median(stimulusRate([stimuli_parameters.Stm.Amplitude] == uAmp(amplitude), cluster))';
    end

    % plot
    plot(uAmp, mean(condition(cluster, :) - condition(cluster, 1), 1), ':o');
    hold on

end

xticks(uAmp)
xticklabels(uAmp(1:nAmp));
xlabel('Pressure (V)')
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

%save figure
figname = sprintf('M%.2i_S%.2i_%s_mean', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

% plot summary fig: boxplot, all units
close all
condition = nan(nClusters, nAmp);

% select data
for amplitude = 1:nAmp
    condition(:, amplitude) = median(stimulusRate([stimuli_parameters.Stm.Amplitude] == uAmp(amplitude), :))';
end

% plot
figure
boxplot(condition)

xticklabels(uAmp(1:nAmp));
xlabel('Pressure (V)')
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

%save figure
figname = sprintf('M%.2i_S%.2i_%s_boxplot_all units', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

%% plot vibrotac stimulus signal
% include also feedback signal
% see plotting in GUI
% not correct yet, Ramp and offset

Waveform = 'BiSine';
Fs = 30000;
StimDur = 500;
Amplitude = 0.5;
SomFreq = 10;
Ramp = 1;
ISI = 1000;
Offset = 0.1;

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

% padding pre-post stimulus time with zeroes
PrePostDur = 0.2 * ISI * 0.001;
PrePostSamp = round(PrePostDur * Fs);
tt = [-flip(1:PrePostSamp)./Fs, tt, tt(end) + (1:PrePostSamp)./Fs];
som_waveform = [zeros(1,PrePostSamp), som_waveform, zeros(1,PrePostSamp)];

% plotting the waveform
figure;
plot(tt(1:length(som_waveform)),som_waveform+Offset);
xlim([min(tt),max(tt)]);

%% -------------------- Local functions -------------------------- %%

% Extract spikes M1
function [spiketimes, cids, cpos, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, TTLPath, messagesPath, relevant_sessions, rec_samples, Fs)
% Kilosort: post-curation unit extraction
% INPUT - paths to sorted data (cluster_info, table), spike times (vector)
% and matched unit ids (vector), recording time stamps (vector)
% OUTPUT - spike times of each single unit
% based on postcuration in PostCuration_ABW.m

% check if paths contain needed files
if ~isfile([KSPath,'cluster_info.tsv']) || ~isfile([KSPath,'spike_clusters.npy']) || ~isfile([KSPath,'spike_times.npy'])
    error('Files not found in KSPath.');
elseif ~isfile([TTLPath,'sample_numbers.npy']) ||  ~isfile([TTLPath,'states.npy'])
    error('Files not found in TTLPath.')
end

cluster_info = readtable([KSPath,'cluster_info.tsv'],'FileType','text'); % info on clusters
cids = cluster_info.cluster_id(strcmp(cluster_info.group,'good'))';
[~, idx] = ismember(cids, cluster_info.cluster_id);
% to do: improve name of table cpos
cpos(:,1) = cluster_info.cluster_id(idx);
cpos(:,2) = cluster_info.ch(idx);
cpos(:,3) = cluster_info.depth(idx);
cpos(:,4) = cluster_info.fr(idx);
cpos(:,5) = cluster_info.n_spikes(idx);
cpos = array2table(cpos, 'VariableNames', {'id' , 'channel', 'depth', 'firing_rate', 'nr_spikes'});

fprintf('Found %i good units for analysis\n', length(cids));

% load files
spike_times = readNPY([KSPath 'spike_times.npy']); % spike_times contains spike time indexing, not time/samplenr itself
spike_clusters = readNPY([KSPath 'spike_clusters.npy']); % matched cluster ids
TTL_samples = readNPY([TTLPath 'sample_numbers.npy']); % sample nr all recorded TTLs
TTL_states = readNPY([TTLPath 'states.npy']);
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

% define and initiate variables
%Fs = 30000; % Sampling frequency (in Hz)
Srise = [];
Sfall = [];

% remove camera TTLs
index = (abs(TTL_states) == 8);
TTL_states(index) = [];
TTL_samples(index) = [];

% recording sessions table
%sessions_TTLs = makeSessionTTLs(messagesPath);
%message_samples = readNPY([messagesPath 'sample_numbers.npy']); % session TTLs
%sessions_TTLs = [];
%sessions_TTLs(:,1) = [1 2 3 3 4 4 5 5 6 6 7 7 8 8 10 10 11 11 12 12 13 13 14 14 17 17]; % session nr
%sessions_TTLs(:,2) = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]; % session start/end (1/0)
%sessions_TTLs(:,3) = message_samples(1:length(sessions_TTLs)); % sample nr
%sessions_TTLs(21:22,:) = []; % remove session 13, aborted halfway

% specify TTL nr to be used
Nr_sessions = (relevant_sessions(1):relevant_sessions(2))';
for file = 1:length(Nr_sessions)
    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);

    Nr_sessions(file,1) = str2double(stimuli_parameters.Par.Set);
    Nr_sessions(file,2) = size(stimuli_parameters.Stm, 1);
    %Nr_sessions(file,2) = str2double(stimuli_parameters.Par.Ntrl);

    % notate TTL nr (TTLS: 5 = SOM, 2 = AUD, 8 = CAM)
    if strcmp(stimuli_parameters.Par.Rec, 'FRA')
        Nr_sessions(file,3) = 2;
    elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
        Nr_sessions(file,3) = 2;
    elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
        Nr_sessions(file,3) = 5;
    end

end

%keep only TTLs recorded during specific session
for i = 1:length(sessions_TTLs)

    disp(['iteration ' num2str(i)])

    if sessions_TTLs(i,2) == 1
        session_start = sessions_TTLs(i,3); %start
        session_end = sessions_TTLs(i+1,3); %end
        idx = (TTL_samples >= session_start) & (TTL_samples < session_end);
        tTTL_states = TTL_states(idx);
        tTTL_samples = TTL_samples(idx);
    elseif (i == 1) && (sessions_TTLs(i,2) == 0) %first session start not noted
        session_start = rec_samples(1); %start
        session_end = sessions_TTLs(i,3); %end
        idx = (TTL_samples >= session_start) & (TTL_samples < session_end);
        tTTL_states = TTL_states(idx);
        tTTL_samples = TTL_samples(idx);
    elseif (i == 2) && (sessions_TTLs(i,2) == 0) % special case in M7
        session_start = sessions_TTLs(i-1,3); %start
        session_end = sessions_TTLs(i,3); %end
        idx = (TTL_samples >= session_start) & (TTL_samples < session_end);
        tTTL_states = TTL_states(idx);
        tTTL_samples = TTL_samples(idx);
    else
        continue
    end

    disp(['session: ' num2str(sessions_TTLs(i,1))])
    disp(['start sample: ' num2str(session_start)])
    disp(['end sample: ' num2str(session_end)])

    % keep only session specific TTLs
    TTLnr = Nr_sessions(sessions_TTLs(i,1), 3); % get TTLnr corresponding to session
    index = (abs(tTTL_states) == TTLnr);
    ttTTL_states = tTTL_states(index);
    ttTTL_samples = tTTL_samples(index);

    % session always starts with high/positive TTL number
    if ttTTL_states(1) < 0
        ttTTL_samples(1) = [];
        ttTTL_states(1) = [];
    end

    % from session keep only the relevant TTLs (so aud ttl for aud session)
    tSrise = ttTTL_samples(ttTTL_states == TTLnr);
    tSfall = ttTTL_samples(ttTTL_states == -TTLnr);

    % remove artefacts where Srise == Sfall
    minDur = 30 ; % samples (= 1ms)
    idx = find ((tSfall - tSrise) < minDur);
    tSrise(idx) = [];
    tSfall(idx) = [];

    disp(['saved TTLs: ' num2str(size(tSrise,1))])

    % output variables
    Srise = [Srise; tSrise];
    Sfall = [Sfall; tSfall];
    %TTLs_sample = [TTLs_sample, tTTL_samples];
    %TTLs_state = [TTLs_state, tTTL_state];

end

% get all spiketimes from each good unit
% to do: make local function
spiketimes = cell(length(cids), 1);

% % exclusion criteria: firingrate > 0.1Hz
total_rec = (rec_samples(length(rec_samples)) - rec_samples(1))/Fs;
minimal_freq = total_rec * 0.1; % min amount of spike

for cluster = 1:length(cids)
    idx = (spike_clusters == cids(cluster));
    spiketimes{cluster} = rec_samples(spike_times(idx));
    % if length(spiketimes{cluster}) < minimal_freq % to do: actually remove unit from all relevant variables
    %     spiketimes{cluster} = NaN;
    % end
end

fprintf('unit extraction done\n');

end



%% CF / threshold / band with 10dB

function FRACfThr(clusterIdx,FRA)

deltaInt = [0,10,20,30,40];
Color =[0,      0,      1; ...
    0.25,   0.25,   1;...
    0.5,    0.5,    1; ...
    0.5,    0.5,    0.5;...
    0,      0,      0] ;%;0.25,0.25,.75;0,0,0];

    for cluster = clusterIdx
        % set alpha value
        alpha = 1e-3;

        % assign significance (based on bootstrapped p-value)
        sig = double(FRA.FACApval(:,:,Set,cluster) <= alpha);
        sig(isnan(FRA.FACApval(:,:,Set,cluster))) = nan;

        % find minimum threshold and CF
        [row,col] = find(sig==1);
        if isempty(row)
            minThr = NaN;
        else
            minThr = FRA.UInt(min(row));
        end
        CF = mean(FRA.UFreq(col(row == min(row))));

        % find BF
        freqMat = repmat(FRA.UFreq',FRA.NInt,1);
        IntMat = repmat(FRA.UInt,1,FRA.NFreq);
        [~,I] = max(FRA.FRASR(:,:,Set,cluster),[],'all','linear');
        if (sig(I) == 1) % condition on significant response
            BF = freqMat(I); BFInt = IntMat(I);
        else
            BF = NaN; BFInt = NaN;
        end


        % Bandwidth
        spontRate = FRA.FRASR(1,1,Set,cluster);
        excBand = NaN(length(deltaInt),2);
        inhBand = NaN(length(deltaInt),2);
        excBW_oct = NaN(length(deltaInt),1);
        inhBW_oct = NaN(length(deltaInt),1);
        excQ = NaN(length(deltaInt),1);
        inhQ = NaN(length(deltaInt),1);
        octStep = 4;
        for k = 1:length(deltaInt)
            i = find(FRA.UInt == minThr+deltaInt(k));
            excBandIdx = (FRA.FRASR(i,:,Set,cluster) >= spontRate) & ...
                (sig(i,:) == 1 );
            inhBandIdx = (FRA.FRASR(i,:,Set,cluster) <= spontRate) & ...
                (sig(i,:) == 1 );

            if(sum(excBandIdx)>0)
                %             excBW(k) = (find(excBandIdx,1,'last') - find(excBandIdx,1,'first') +1) / octStep;
                xx = find(excBandIdx,1,'first');
                excBand(k,1) = sqrt(FRA.UFreq(xx) * FRA.UFreq(xx-1)) ;
                xx = find(excBandIdx,1,'last');
                excBand(k,2) = sqrt(FRA.UFreq(xx) * FRA.UFreq(xx+1)) ;

                excBW_oct(k) = log2(excBand(k,2) / excBand(k,1));
                excQ(k) = CF / (excBand(k,2) - excBand(k,1) );
            end

            if(sum(inhBandIdx)>0)
                %             inhBW(k) = (find(inhBandIdx,1,'last') - find(inhBandIdx,1,'first') +1) / octStep;
                xx = find(inhBandIdx,1,'first');
                inhBand(k,1) = sqrt(FRA.UFreq(xx) * FRA.UFreq(xx-1)) ;
                xx = find(inhBandIdx,1,'last');
                inhBand(k,2) = sqrt(FRA.UFreq(xx) * FRA.UFreq(xx+1)) ;

                inhBW_oct(k) = log2(inhBand(k,2) / inhBand(k,1));
                inhQ(k) = CF / (inhBand(k,2) - inhBand(k,1) );
            end
        end
        factor = FRA.FRASR(1,1,Set,cluster)/FRA.FRAScnt(1,1,Set,cluster);
        spontRateSD = factor*FRA.FRAScntSD(1,1,Set,cluster);


        % evaluate output
        fig = figure(setNum*100+cluster);clf;
        set(fig,'Position',[100,100,1600,800]);

        % spike rate and p-values
        h = subplot(2,3,1);

        CData = FRA.FRASR(:,:,Set,cluster);
        imagesc(h,CData,'AlphaData',~isnan(CData),[0,Inf]);
        % contour of FACA p-value
        Cont = -log10(FRA.FACApval(:,:,Set,cluster));
        hold(h,'on');contour(h,Cont,[2,3],'w','ShowText','on');hold(h,'off');
        % contour of max neighbour correlation
        Cont = (FRA.MaxNeighCorr(:,:,Set,cluster));
        hold(h,'on');contour(h,Cont,[0.1,0.2,0.5],'r','ShowText','on');hold(h,'off');
        % format and label graph
        set(h,'Xscale','lin','YDir','normal',...
            'FontName','Arial','FontWeight','bold','FontSize',12, ...
            'XTick',2:4:FRA.NFreq,'XTickLabel',round(FRA.UFreq(2:4:FRA.NFreq),1), 'XTickLabelRotation',45,...
            'YTick',2:2:FRA.NInt,'YTicklabel',FRA.UInt(2:2:FRA.NInt));
        ylabel(h,'Intensity (dB SPL)');
        xlabel(h,'frequency (kHz)')
        title(['BF = ',num2str(BF,'%.1f'),' kHz  ',' Int = ',num2str(BFInt),' dB'])
        cb = colorbar(h,'eastoutside');
        cb.Label.String = 'spike rate (spk/s)';

        % thresholded significant response
        h = subplot(2,3,2);

        imagesc(h, sig,'AlphaData',~isnan(sig),[0,1]);
        set(h,'Xscale','lin','YDir','normal',...
            'FontName','Arial','FontWeight','bold','FontSize',12, ...
            'XTick',2:4:FRA.NFreq,'XTickLabel',round(FRA.UFreq(2:4:FRA.NFreq),1), 'XTickLabelRotation',45,...
            'YTick',2:2:FRA.NInt,'YTicklabel',FRA.UInt(2:2:FRA.NInt));
        title(['CF = ',num2str(CF,'%.1f'),' kHz  ',' minThr = ',num2str(minThr),' dB'])
        ylabel(h,'Intensity (dB SPL)');
        xlabel(h,'frequency (kHz)')
        cb = colorbar(h,'eastoutside','Ticks',[0,1],'TickLabels',{'not sig','sig'});
        cb.Label.String = 'significance';

        % MED FSL
        h = subplot(2,3,3);

        imagesc(h, FRA.MedFSL(:,:,Set,cluster)*1e3,'AlphaData',~isnan(sig),[5,40]);
        set(h,'Xscale','lin','YDir','normal',...
            'FontName','Arial','FontWeight','bold','FontSize',12, ...
            'XTick',2:4:FRA.NFreq,'XTickLabel',round(FRA.UFreq(2:4:FRA.NFreq),1), 'XTickLabelRotation',45,...
            'YTick',2:2:FRA.NInt,'YTicklabel',FRA.UInt(2:2:FRA.NInt));
        BF_FSL = FRA.MedFSL(FRA.UInt==BFInt,FRA.UFreq==BF,Set,cluster)*1e3;
        [~,CFIdx] = min(abs(FRA.UFreq-CF));
        CF30_FSL = FRA.MedFSL(FRA.UInt==minThr+30,CFIdx,Set,cluster)*1e3;
        minFSL = min(FRA.MedFSL(:,:,Set,cluster)*1e3,[],'all','omitnan');
        title(['@CF (Thr+30dB): ',num2str(CF30_FSL,'%.1f'),' ms  ',newline,' @BF: ',num2str(BF_FSL,'%.1f'),' ms  ' ...
            ,newline,' min: ',num2str(minFSL,'%.1f'), ' ms'])
        cb = colorbar(h,'eastoutside');
        cb.Label.String = 'median first spike latency (ms)';

        ylabel(h,'Intensity (dB SPL)');
        xlabel(h,'frequency (kHz)')

        % bandwidth
        h = subplot(2,3,4);cla(h);hold(h,'on');
        lbls = {};traces =[];
        for k = 1:length(deltaInt)
            i = find(FRA.UInt == minThr+deltaInt(k));
            if deltaInt(k) == 0; lbl = 'Thr'; else; lbl = ['Thr + ' num2str(deltaInt(k)) ' dB'];end
            lbls = [lbls,{lbl}];
            traces(k) = plot(h,FRA.UFreq,FRA.FRASR(i,:,Set,cluster),'Color',Color(k,:));
            sigFreq = sig(i,:) == 1;
            plot(h,FRA.UFreq(sigFreq),FRA.FRASR(i,sigFreq,Set,cluster),'o','Color',Color(k,:),'MarkerFaceColor',Color(k,:));
        end
        %     meanRate = mean(FRA.FRASR(:,:,setIdx,cluster),'all','omitnan');
        traces(k+1) = plot(h,FRA.UFreq,repmat(spontRate,FRA.NFreq,1),'--k');
        lbls = [lbls,{'spont rate'}];
        %     traces(k+2) = plot(h,FRA.UFreq,repmat(spontRate+3*spontRateSD/sqrt(20),FRA.NFreq,1),':b');
        %     lbls = [lbls,{'spont + 3 sem'}];
        set(h,'Xscale','log',...
            'FontName','Arial','FontWeight','bold','FontSize',12, ...
            'XTick',FRA.UFreq(2:4:FRA.NFreq),'XTickLabel',round(FRA.UFreq(2:4:FRA.NFreq),1), 'XTickLabelRotation',45 ...
            );
        ylabel(h,'spike rate (spk/s)');ylim([0,inf]);
        xlabel(h,'frequency (kHz)')
        legend(traces,lbls,'location','eastoutside')
        title(['alpha = ',num2str(alpha,'%.4f'),''])

        % bandwidth 2
        h = subplot(2,3,5);cla(h);hold(h,'on');
        traces = [];
        traces(2) = plot(deltaInt+1,inhBW_oct,'-v','Color',[.75,0,.75],'MarkerFaceColor',[.75,0,.75]);
        traces(1) = plot(deltaInt-1,excBW_oct,'-square','Color',[0,.5,0],'MarkerFaceColor',[0,.5,0]);
        set(h,'FontName','Arial','FontWeight','bold','FontSize',12 ...
            ,'XTick',deltaInt,'XTickLabel',deltaInt);
        ylabel('bandwidth (octave)');ylim([0,2.5]);
        xlabel('Intensity (dB re. minThr)');xlim([min(deltaInt)-2.5,max(deltaInt)+2.5]);
        legend(traces,{'excitatory','inhibitory'},'Location','eastoutside');

        % bandwidth 3
        h = subplot(2,3,6);cla(h);hold(h,'on');
        traces = [];
        traces(2) = plot(deltaInt+1,inhQ,'-v','Color',[.75,0,.75],'MarkerFaceColor',[.75,0,.75]);
        traces(1) = plot(deltaInt-1,excQ,'-square','Color',[0,.5,0],'MarkerFaceColor',[0,.5,0]);
        set(h,'FontName','Arial','FontWeight','bold','FontSize',12 ...
            ,'XTick',deltaInt,'XTickLabel',deltaInt);
        ylabel('QXdB [CF/bw at X dB re. Thr]');ylim([0,10]);
        xlabel('Intensity (dB re. minThr)');xlim([min(deltaInt)-2.5,max(deltaInt)+2.5]);
        legend(traces,{'excitatory','inhibitory'},'Location','eastoutside');


        sgtitle([Mouse ' - E' num2str(cpos(cluster)) ' - FRA ',num2str(Set),' Set #',num2str(setNum)]);
        figName = [Mouse '- E' num2str(cpos(cluster)) '_FRA',num2str(setNum)];
        set(fig,'Name',figName);

        %     saveas(fig,[ResPath '\Figures\' figName,'.png']);
    end
end
