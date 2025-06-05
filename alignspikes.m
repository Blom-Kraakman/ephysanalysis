function [aligned_spikes] = alignspikes(BehaviorPath, TTLPath, OutPath, spiketimes, relevant_sessions, skip_sessions, cids, sessions_TTLs, Fs)
% align spikes
% INPUT - spiketimes (cell array, 1 cell: timestamp of spike for 1 unit),
% TTLs (vector), stimulus parameters (struct)
% OUTPUT - cell aray (aligned_spikes) where 1 cell contains all spike times of one unit
% relative to stim onset
% based on AlignSpkS in PostCuration_ABW.m

% check input
if isempty(spiketimes) || isempty(cids)
    error('Input arguments missing. Extract spikes first.')
elseif isempty(relevant_sessions)
    relevant_sessions = input('Enter first and last session number in this recording: ');
end

%stim_counter = 0;
aligned_spikes = {}; % nStim x Ncids

% select correct behaviour file
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

% load TTL data
TTL_samples = readNPY([TTLPath 'sample_numbers.npy']); % sample nr all recorded TTLs
TTL_states = readNPY([TTLPath 'states.npy']);
TTL_words = readNPY([TTLPath 'full_words.npy']);

% remove camera TTLs
index = (abs(TTL_states) == 8);
TTL_states(index) = [];
TTL_samples(index) = [];
TTL_words(index) = [];


% allign spikes to stimuli (per session)
for file = relevant_sessions(1):relevant_sessions(2)

    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);
   
    % check file selection
    if file ~= str2double(stimuli_parameters.Par.Set)
        warning('WARNING: File indexing error, check if correct stimulus file is being selected.')
        break
    end

    % skip earlier stim files not in rec session
    if ismember(str2double(stimuli_parameters.Par.Set), skip_sessions) || ~ismember(str2double(stimuli_parameters.Par.Set), relevant_sessions(1):relevant_sessions(2))
        continue
    else
        disp(['Aligning session: ' num2str(file)])
    end

    % select correct analysis window and
    if strcmp(stimuli_parameters.Par.Rec, 'SOM')
        PreT = str2double(stimuli_parameters.Par.SomatosensoryISI)/4; % amount of msec. to include before Srise;
        PostT = (str2double(stimuli_parameters.Par.SomatosensoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI)/4); % amount of msec. to include after Sfall;
    % -- NEW -- %
    elseif strcmp(stimuli_parameters.Par.Rec, 'SxA') % TO TEST
        if length(str2num(stimuli_parameters.SomAudSOA)) > 2 % multiple SOA delays
            PreT = str2double(stimuli_parameters.Par.SomatosensoryISI)/2;
            % PostT = (str2double(stimuli_parameters.Par.AuditoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI)/2);
            PostT = (min(str2double(stimuli_parameters.Par.AuditoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryStimTime) + max(stimuli_parameters.Stm.SomAudSOA) ...
                , str2double(stimuli_parameters.Par.SomatosensoryISI)/2));
        else
    % -- % 
            PreT = str2double(stimuli_parameters.Par.SomatosensoryISI)/2;
            % PostT = (str2double(stimuli_parameters.Par.AuditoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI)/2);
            PostT = (max(str2double(stimuli_parameters.Par.AuditoryStimTime), str2double(stimuli_parameters.Par.SomatosensoryStimTime)) ...
                + str2double(stimuli_parameters.Par.SomatosensoryISI)/2); % take max stim time
        end
    elseif strcmp(stimuli_parameters.Par.Rec, 'FRA')
        PreT  = str2double(stimuli_parameters.Par.FRAStimTime);
        PostT = str2double(stimuli_parameters.Par.FRAPostTime);
    elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
        PreT = (str2double(stimuli_parameters.Par.AMPreTime) + (str2double(stimuli_parameters.Par.AMPostTime)/4));
        PostT = (str2double(stimuli_parameters.Par.AMStimTime) + (str2double(stimuli_parameters.Par.AMPostTime)/2));
    end

    % select correct number of TTLs based on stimuli_parameters.par file
    SpkT = [];
    NStim = size(stimuli_parameters.Stm, 1); % nr trials to align per stim session
    disp(['Stimuli in session: ' num2str(NStim)])

    % get Srise & Sfall for session
    tsessions_TTLs = sessions_TTLs(sessions_TTLs(:,1) == file, :);
    [Srise, Sfall] = getTrialTTLs(tsessions_TTLs, TTL_states, TTL_samples, TTL_words);
    Srise = double(Srise);
    Sfall = double(Sfall);

    if NStim ~= size(Srise)
        warning('Check Srise')
        continue
    end

    % loop through all clusters
    for cluster = 1:length(cids)

        % pre-allocate output cell array
        tSpkT = cell(NStim, 1);
        Spks = double(spiketimes{cluster});

        % loop through all stimuli
        for stim = 1:NStim
            %tS = (Spks - Srise(stim + stim_counter))./ Fs; % s; spiketimes re. stimulus k onset
            tS = (Spks - Srise(stim))./ Fs; % s; all spiketimes re. stimulus k onset
            sel = (tS > - PreT/1e3) & (tS < PostT/1e3); % select relevant spikes (including pre and post rec)
            tSpkT{(stim)} = tS(sel); % store spike times in cell array
        end

        SpkT = [SpkT, tSpkT]; % stimuli x units per session

    end

    %stim_counter = stim_counter + NStim; % keep track of how many stimuliu have past
    aligned_spikes = [aligned_spikes; SpkT]; % stimuli x units recording

    % save aligned spikes
   filename = sprintf('M%.2i_S%.2i_%s_AlignedSpikes', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
   save(fullfile(OutPath, filename), "SpkT","Srise","Sfall")

end

fprintf('spike alignment done\n');

end
