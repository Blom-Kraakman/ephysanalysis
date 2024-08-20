function [aligned_spikes] = alignspikes(BehaviorPath, OutPath, spiketimes, relevant_sessions, skip_sessions, Srise, Sfall, cids, Fs)
% align spikes
% INPUT - spiketimes (cell array, 1 cell: timestamp of spike for 1 unit),
% TTLs (vector), stimulus parameters (struct)
% OUTPUT - cell aray (aligned_spikes) where 1 cell contains all spike times of one unit
% relative to stim onset
% based on AlignSpkS in PostCuration_ABW.m

% check input
if isempty(spiketimes) || isempty(cids) || isempty(Srise) || isempty(Sfall)
    error('Input arguments missing. Extract spikes first.')
elseif isempty(relevant_sessions)
    relevant_sessions = input('Enter first and last session number in this recording: ');
end

Srise = double(Srise);
Sfall = double(Sfall);
stim_counter = 0;
aligned_spikes = {}; % nStim x Ncids

% select correct behaviour file
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

% allign spikes to stimuli (per session)
for file = relevant_sessions(1):relevant_sessions(2)

    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);

    if ismember(str2double(stimuli_parameters.Par.Set), skip_sessions)
        continue
    end

    % select correct analysis window and
    if strcmp(stimuli_parameters.Par.Rec, 'SOM')
        PreT = str2double(stimuli_parameters.Par.SomatosensoryISI)/4; % amount of msec. to include before Srise;
        PostT = (str2double(stimuli_parameters.Par.SomatosensoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI)/4); % amount of msec. to include after Sfall;
    elseif strcmp(stimuli_parameters.Par.Rec, 'SxA')
        PreT = str2double(stimuli_parameters.Par.SomatosensoryISI)/2;
        %PostT = (str2double(stimuli_parameters.Par.AuditoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI)/2);
        PostT = (max(str2double(stimuli_parameters.Par.AuditoryStimTime), str2double(stimuli_parameters.Par.SomatosensoryStimTime)) ...
            + str2double(stimuli_parameters.Par.SomatosensoryISI)/2); % take max stim time
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

    disp(['session ' num2str(stimuli_parameters.Par.Set)])
    disp(['Nr stim in session ' num2str(NStim)])

    % loop through all clusters
    for cluster = 1:length(cids)

        % pre-allocate output cell array
        tSpkT = cell(NStim, 1);
        Spks = double(spiketimes{cluster});

        % loop through all stimuli
        for stim = 1:NStim
            tS = (Spks - Srise(stim + stim_counter))./ Fs; % s; spiketimes re. stimulus k onset
            sel = (tS > - PreT/1e3) & (tS < PostT/1e3); % select relevant spikes (including pre and post rec)
            tSpkT{(stim)} = tS(sel); % store spike times in cell array
        end

        % if strcmp(stimuli_parameters.Par.Rec, 'SxA')
        %     % add delay to "SO","OO"
        %     idx_Sound = find(ismember(stimuli_parameters.Stm.MMType,["SA","OA"]));
        %     idx_noSound = find(ismember(stimuli_parameters.Stm.MMType,["SO","OO"]));
        %     for ii = idx_noSound'
        %         aligned_spikes{ii,cluster} = aligned_spikes{ii,cluster} - 0.25;
        %     end
        % end

        SpkT = [SpkT, tSpkT]; % stimuli x units per session

    end

    stim_counter = stim_counter + NStim; % keep track of how many stimuliu have past
    aligned_spikes = [aligned_spikes; SpkT]; % stimuli x units recording

    disp(['stim counter: ' num2str(stim_counter)])

    % save aligned spikes
    filename = sprintf('M%.2i_S%.2i_%s_AlignedSpikes', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
    save(fullfile(OutPath, filename), "SpkT","Srise","Sfall")

end

fprintf('spike alignment done\n');

end
