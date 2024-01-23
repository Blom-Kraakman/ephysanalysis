function [aligned_spikes] = alignspikes(BehaviorPath, spiketimes, relevant_sessions, Srise, Sfall, cids)
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

Fs = 30000; % Sampling frequency (in Hz)
Srise = double(Srise);
Sfall = double(Sfall);
stim_counter = 0;
aligned_spikes = {}; % nStim x Ncids

% select correct behaviour file
stim_files = dir(fullfile(BehaviorPath, '\*.mat'));

% allign spikes to stimuli (per session)
for file = relevant_sessions(1):relevant_sessions(2)

    stimuli_parameters = load([stim_files(file).folder '\' stim_files(file).name]);

    % select correct analysis window and
    if strcmp(stimuli_parameters.Par.Rec, 'SOM')
        PreT  = 200; % amount of msec. to include before Srise;
        PostT = 700; % amount of msec. to include after Sfall;
    elseif strcmp(stimuli_parameters.Par.Rec, 'FRA')
        PreT  = str2double(stimuli_parameters.Par.FRAStimTime);
        PostT = str2double(stimuli_parameters.Par.FRAPostTime);
    elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
        PreT = str2double(stimuli_parameters.Par.AMStimTime);
        PostT = str2double(stimuli_parameters.Par.AMPostTime);
        % AMStimTime 200, AMPostTime 400, AMPreTime 100
    end

    % select correct number of TTLs based on stimuli_parameters.par file
    NStim = size(stimuli_parameters.Stm, 1); % nr trials to align per stim session
    SpkT = [];

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

        SpkT = [SpkT, tSpkT]; % stimuli x units per session

    end

    stim_counter = stim_counter + NStim; % keep track of how many stimuliu have past
    aligned_spikes = [aligned_spikes; SpkT]; % stimuli x units recording

    % save aligned spikes
    set = sprintf('%02d', str2double(stimuli_parameters.Par.Set));
    filename = ['M04_S' set '_' stimuli_parameters.Par.Rec '_AlignedSpikes'];
    save(fullfile('D:\DATA\Processed', filename), "SpkT")

end

fprintf('spike alignment done\n');

end
