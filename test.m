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



%%
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

%%
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


%% Extract spikes M1
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

%% plot SOM piezo data
% if strcmp(stimuli_parameters.Par.SomatosensoryActuator, 'Piezo')
%     %if strcmp(stimuli_parameters.Par.SomatosensoryActuator, 'Piezo') || strcmp(stimuli_parameters.Par.Rec, 'SOM')
% 
%     xrange = [-.2, +.7];
% 
%     for cluster = 1:length(cids)
% 
%         figure;
% 
%         fig = subplot(2,1,1); % rasterplot
%         set(gcf,'position',[500,150,800,600])
%         %set(gcf,'position',[500,150,600,400])
% 
%         % make rasterplot
%         Var = [stimuli_parameters.Stm.Amplitude, stimuli_parameters.Stm.SomFreq];
%         %[f, YTick, YTickLab] = plotraster(fig, aligned_spikes(1:40, 1), Var, [0, 0, 0], [], 1);
% 
%         [f, YTick, YTickLab, varargout] = plotraster(fig, aligned_spikes(:, cluster), Var, [0, 0, 0], [], 1);
%         yticks(YTick{1});
%         yrange = [min(YTick{2}) - 5, max(YTick{2}) + 5];
%         ylim(f,yrange);
%         yticklabels(unique(stimuli_parameters.Stm.Amplitude));
%         xlim(f,xrange);
% 
%         % format axis
%         xlabel('Time (s)')
%         ylabel('Stimulus off / on')
%         fig.FontSize = 11;
%         title(['Cluster ' num2str(cids(cluster)) ' - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])
% 
%         % make histogram / PSTH
%         SOM_idx = (stimuli_parameters.Stm.Amplitude == 10); % select stimulus trials
%         ctrl_idx = (stimuli_parameters.Stm.Amplitude == 0); % select control trials
%         preT  = -0.2;
%         postT = 0.7;
%         binsize = 0.5;
% 
%         fig = subplot(2,1,2);
% 
%         [N, edges] = histcounts(vertcat(aligned_spikes{SOM_idx, cluster}), preT:binsize:postT);
%         % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 1))/binsize), 'FaceColor', '#D95319') % spike/s
%         plot(edges(1:end-1),((N/sum(SOM_idx == 1))/binsize),'Color', '#D95319','LineWidth',2.5)
%         hold on
%         [N,edges] = histcounts(vertcat(aligned_spikes{ctrl_idx, cluster}), preT:binsize:postT);
%         % histogram('BinEdges', edges, 'BinCounts', ((N/sum(SOM_idx == 0))/binsize), 'FaceColor', '#0072BD')
%         plot(edges(1:end-1),((N/sum(SOM_idx == 0))/binsize),'Color', '#0072BD','LineWidth',2.5)
% 
%         %format axis
%         legend('stimulus', 'control')
%         xlabel('Time (s)')
%         ylabel('Spike rate (Hz)')
%         fig.FontSize = 11;
%         xlim(fig,xrange);
% 
%         sgtitle(['Cluster ' num2str(cids(cluster))])
% 
%         % save plot
%         figname = sprintf('M%.2i_S%.2i_%s_cluster_%i', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec, cids(cluster));
%         saveas(gcf, fullfile(OutPath, [figname '.jpg']));
%         saveas(fig, fullfile(OutPath, figname));
% 
%         hold off
% 
%     end
% end