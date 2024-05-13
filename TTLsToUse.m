function [Srise, Sfall] = TTLsToUse(sessions_TTLs, TTLPath, rec_samples)

% load data
TTL_samples = readNPY([TTLPath 'sample_numbers.npy']); % sample nr all recorded TTLs
TTL_states = readNPY([TTLPath 'states.npy']);
TTL_words = readNPY([TTLPath 'full_words.npy']);

% remove camera TTLs
index = (abs(TTL_states) == 8);
TTL_states(index) = [];
TTL_samples(index) = [];
TTL_words(index) = [];

% define and initiate variables
Srise = [];
Sfall = [];

for i = 1:size(sessions_TTLs, 1)

    disp(['iteration ' num2str(i)])

    % define start and end of session
    if sessions_TTLs(i,2) == 1
        session_start = sessions_TTLs(i,3); %start
        session_end = sessions_TTLs(i+1,3); %end
    elseif (i == 1) && (sessions_TTLs(i,2) == 0) % missing start first session
        session_start = rec_samples(1); %start
        session_end = sessions_TTLs(i,3); %end
    elseif (i == 9) && (sessions_TTLs(i,2) == 0) % special case in M6
        session_start = sessions_TTLs(i,3); %start
        session_end = rec_samples(end);
        % elseif (i == 2) && (sessions_TTLs(i,2) == 0) % special case in M7
        %         session_start = sessions_TTLs(i-1,3); %start
        %         session_end = sessions_TTLs(i,3); %end
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

    %special case M7
    % if sessions_TTLs(i, 1) == 14
    %     idx = (abs(tTTL_states) == 2);
    %     tTTL_states = tTTL_states(idx);
    %     tTTL_words = tTTL_words(idx);
    %     tTTL_samples = tTTL_samples(idx);
    %     clear idx;
    % end

    % keep only session specific TTLs, gives 1 when either and both stim (2: aud, 16(2^4): som) are on/high
    stim_on = bitand(tTTL_words, (2+16)) > 0; % numbers depend on session type
    stim_rise = [true; diff(stim_on) > 0];
    stim_fall = [false; diff(stim_on) < 0];

    % group TTLs on rising or falling edge
    tSrise = tTTL_samples(stim_rise);
    tSfall = tTTL_samples(stim_fall);

    idx = (tSrise(1) == tTTL_samples(3));
    tSrise(idx) = [];
    
    % remove artefacts where Srise == Sfall
    minDur = 30; % samples (= 1ms)
    idx = find((tSfall - tSrise) < minDur);
    tSrise(idx) = [];
    tSfall(idx) = [];

    disp(['saved TTLs: ' (num2str(size(tSrise,1)*2))])

    % output variables
    Srise = [Srise; tSrise];
    Sfall = [Sfall; tSfall];

end
end
