function [spiketimes, cids, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, TTLPath, messagesPath, relevant_sessions, rec_samples, Fs, OutPath)
% Kilosort: post-curation unit extraction
% INPUT - paths to sorted data (cluster_info, table), spike times (vector)
% and matched unit ids (vector), recording time stamps (vector)
% OUTPUT - spike times of each single unit
% based on postcuration in PostCuration_ABW.m

%check if paths contain needed files
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
TTL_words = readNPY([TTLPath 'full_words.npy']);
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

% remove camera TTLs
index = (abs(TTL_states) == 8);
TTL_states(index) = [];
TTL_samples(index) = [];
TTL_words(index) = [];

% recording sessions table
sessions_TTLs = getSessionTTLs(messagesPath, rec_samples, Fs);

% specify TTL nr to be used
Nr_sessions = (relevant_sessions(1):relevant_sessions(2))';
for file = 1:length(Nr_sessions)
    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);

    Nr_sessions(file,1) = str2double(stimuli_parameters.Par.Set);
    Nr_sessions(file,2) = stimuli_parameters.Par.Ntrl;

    % notate TTL nr (TTLS: 5 = SOM, 2 = AUD, 8 = CAM)
    if strcmp(stimuli_parameters.Par.Rec, 'FRA')
        Nr_sessions(file,3) = 2;
    elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
        Nr_sessions(file,3) = 2;
    elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
        Nr_sessions(file,3) = 5;
    elseif strcmp(stimuli_parameters.Par.Rec, 'SxA')
        Nr_sessions(file,3) = 0;
    else
        error('Did not recognize session type.')
    end

end

% remove unexecuted sessions
sessions_TTLs([1,9,21],:) = [];
Nr_sessions(10,:) = [];

%keep only TTLs recorded during specific session
[Srise, Sfall] = TTLsToUse(sessions_TTLs, TTL_samples, TTL_states, TTL_words, rec_samples, Nr_sessions);

% get all spiketimes from each good unit
spiketimes = cell(length(cids), 1);

% exclusion criteria: firingrate > 0.1Hz
total_rec = (rec_samples(length(rec_samples)) - rec_samples(1))/Fs;
minimal_freq = total_rec * 0.1; % min amount of spike

for cluster = 1:length(cids)
    idx = (spike_clusters == cids(cluster));
    spiketimes{cluster} = rec_samples(spike_times(idx));
    % if length(spiketimes{cluster}) < minimal_freq % to do: actually remove unit from all relevant variables
    %     spiketimes{cluster} = NaN;
    % end
end

filename = sprintf('M%.2i_S%02d-%02d_InfoGoodUnits', str2double(stimuli_parameters.Par.MouseNum), relevant_sessions(1), relevant_sessions(2));
save(fullfile(OutPath, filename), "cpos") %cpos variables: unit id, channel, depth, avg firing rate, nr spikes;

fprintf('unit extraction done\n');

end

function [Srise, Sfall] = TTLsToUse(sessions_TTLs, TTL_samples, TTL_states, TTL_words, rec_samples, Nr_sessions)

% define and initiate variables
Srise = [];
Sfall = [];

for i = 1:length(sessions_TTLs)

    disp(['iteration ' num2str(i)])

    % define start and end of session
    if sessions_TTLs(i,2) == 1
        session_start = sessions_TTLs(i,3); %start
        session_end = sessions_TTLs(i+1,3); %end
    elseif (i == 1) && (sessions_TTLs(i,2) == 0) % missing start first session 
        session_start = rec_samples(1); %start
        session_end = sessions_TTLs(i,3); %end
        % elseif (i == 9) && (sessions_TTLs(i,2) == 0) % special case in M6
        %     session_start = sessions_TTLs(i,3); %start
        %     session_end = rec_samples(end);
        % elseif (i == 2) && (sessions_TTLs(i,2) == 0) % special case in M7
        %     session_start = sessions_TTLs(i-1,3); %start
        %     session_end = sessions_TTLs(i,3); %end
    else
        continue
    end

    % retrieve all session samples
    idx = (TTL_samples >= session_start) & (TTL_samples < session_end);
    tTTL_states = TTL_states(idx);
    tTTL_words = TTL_words(idx);
    tTTL_samples = TTL_samples(idx);
    clear idx;

    disp(['session: ' num2str(sessions_TTLs(i,1))])
    disp(['start sample: ' num2str(session_start)])
    disp(['end sample: ' num2str(session_end)])
    disp(['Session related samples: ' num2str(size(tTTL_states, 1))])

    % SOMETHING GOING WRONG WITH FLICKER AT START
    % keep only session specific TTLs
    % gives 1 when either and both stim (2: aud, 16(2^4): som) are on/high
    stim_on = bitand(tTTL_words, (2+16)) > 0; % numbers depend on session type?
    stim_rise = [true; diff(stim_on) > 0];
    stim_fall = [false; diff(stim_on) < 0];

    % idx = (Nr_sessions(:,1) == sessions_TTLs(i,1)); % get TTLnr corresponding to session
    % TTLnr = Nr_sessions(idx, 3); %TTLnr = Nr_sessions(sessions_TTLs(i,1), 3);
    % index = (abs(tTTL_states) == TTLnr);
    %
    % % SxA session TTLs
    % if Nr_sessions(idx, 3) == 0
    %     TTLnr = [2, 5];
    %     index = ((abs(tTTL_states) == TTLnr(1)) | (abs(tTTL_states) == TTLnr(2)));
    % end

    %ttTTL_states = tTTL_states(index);
    %ttTTL_samples = tTTL_samples(index);

    % % session always starts with high/positive TTL number
    % if ttTTL_states(1) < 0
    %     ttTTL_samples(1) = [];
    %     ttTTL_states(1) = [];
    % end

    % % group TTLs on rising or falling edge
    % tSrise = ttTTL_samples(ttTTL_states >= 0);
    % tSfall = ttTTL_samples(ttTTL_states <= 0);
    
    % group TTLs on rising or falling edge
    tSrise = tTTL_samples(stim_rise);
    tSfall = tTTL_samples(stim_fall);

    idx = (tSrise(1) == tTTL_samples(2));
    tSrise(idx) = [];
    
    % remove artefacts where Srise == Sfall
    minDur = 30; % samples (= 1ms)
    idx = find ((tSfall - tSrise) < minDur);
    tSrise(idx) = [];
    tSfall(idx) = [];

    disp(['saved TTLs: ' (num2str(size(tSrise,1)*2))])

    % output variables
    Srise = [Srise; tSrise];
    Sfall = [Sfall; tSfall];
    %TTLs_sample = [TTLs_sample, tTTL_samples];
    %TTLs_state = [TTLs_state, tTTL_state];

end
end