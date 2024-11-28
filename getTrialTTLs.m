function [Srise, Sfall] = getTrialTTLs(tsessions_TTLs, TTL_states, TTL_samples, TTL_words)

% define start and end of session
if tsessions_TTLs(1,2) == 1 && tsessions_TTLs(2,2) == 0
    session_start = tsessions_TTLs(1,3); %start TTL
    session_end = tsessions_TTLs(2,3); %end TTL
else
    warning('Check sessions_TTLs file')
end

% retrieve all session samples
idx = (TTL_samples >= session_start) & (TTL_samples < session_end);
tTTL_states = TTL_states(idx);
tTTL_words = TTL_words(idx);
tTTL_samples = TTL_samples(idx);
%clear idx;

disp([' start sample: ' num2str(session_start)])
disp([' end sample: ' num2str(session_end)])
disp([' total session related samples: ' num2str(size(tTTL_states, 1))])

% keep only session specific TTLs, gives 1 when either and both stim (2: aud, 16(2^4): som) are on/high
stim_on = bitand(tTTL_words, (2+16)) > 0; % numbers depend on session type
stim_rise = [true; diff(stim_on) > 0];
stim_fall = [false; diff(stim_on) < 0];

% group TTLs on rising or falling edge
Srise = tTTL_samples(stim_rise);
Sfall = tTTL_samples(stim_fall);

% remove artefacts
idx = (Srise(1) == tTTL_samples(3));
Srise(idx) = [];

% remove artefacts where Srise == Sfall
minDur = 30; % samples (= 1ms)
idx = find((Sfall - Srise) < minDur);
Srise(idx) = [];
Sfall(idx) = [];

disp(['Saved TTLs: ' (num2str(size(Srise,1)*2))])

end
