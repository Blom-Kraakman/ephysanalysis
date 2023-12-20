function [spiketimes, cids, cpos, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, recPath, TTLPath, messagesPath, relevant_sessions)
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
elseif ~isfile([recPath, 'sample_numbers.npy'])
    error('Files not found in RecPath.')
end

cluster_info = readtable([KSPath,'cluster_info.tsv'],'FileType','text'); % info on clusters
cids = cluster_info.cluster_id(strcmp(cluster_info.group,'good'))';
[~, idx] = ismember(cids, cluster_info.cluster_id);
cpos(:,1) = cluster_info.cluster_id(idx);
cpos(:,2) = cluster_info.ch(idx);
cpos(:,3) = cluster_info.depth(idx);
cpos(:,4) = cluster_info.fr(idx);
cpos(:,5) = cluster_info.n_spikes(idx);
cpos = array2table(cpos, 'VariableNames', {'id' , 'channel', 'depth', 'firing_rate', 'nr_spikes'});

fprintf('Found %i good units for analysis\n', length(cids));

% load files
rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
spike_times = readNPY([KSPath 'spike_times.npy']); % spike_times contains spike time indexing, not time/samplenr itself
spike_clusters = readNPY([KSPath 'spike_clusters.npy']); % matched cluster ids
message_samples = readNPY([messagesPath 'sample_numbers.npy']); % session TTLs
TTL_samples = readNPY([TTLPath 'sample_numbers.npy']); % sample nr all recorded TTLs
TTL_states = readNPY([TTLPath 'states.npy']);
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

Fs = 30000; % Sampling frequency (in Hz)

% remove camera TTLs
index = (abs(TTL_states) == 8);
TTL_states(index) = [];
TTL_samples(index) = [];

% recording sessions table
% to do: make imported csv
sessions_TTLs(:,1) = [1 1 2 2 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 19 19 20 20 21 21 22 22 23 23 24 24 25 25 26 26 27 27]; % session nr
sessions_TTLs(:,2) = [1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]; % session start/end (1/0)
sessions_TTLs(:,3) = message_samples; % sample nr

% make sessions type + nr trials table
% --> make into function?
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

Srise = [];
Sfall = [];

%keep only TTLs recorded during specific session
for i = 1:length(sessions_TTLs)

    if sessions_TTLs(i,2) == 1
        session_start = sessions_TTLs(i,3); %start
        session_end = sessions_TTLs(i+1,3); %end
        idx = (TTL_samples >= session_start) & (TTL_samples < session_end);
        tTTL_states = TTL_states(idx);
        tTTL_samples = TTL_samples(idx);
    else
        continue
    end

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

    % output variables
    Srise = [Srise, tSrise];
    Sfall = [Sfall, tSfall];
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

fprintf('Unit extraction done\n');

end
