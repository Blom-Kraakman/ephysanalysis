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

