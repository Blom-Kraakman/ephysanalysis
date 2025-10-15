function [SpkT] = alignspikes(spiketimes, cids, Srise, stimuli_parameters, Fs)
% align spikes
% INPUT - spiketimes (cell array, 1 cell: timestamp of spike for 1 unit),
% TTLs (vector), stimulus parameters (struct)
% OUTPUT - cell aray (aligned_spikes) where 1 cell contains all spike times of one unit
% relative to stim onset
% based on AlignSpkS in PostCuration_ABW.m

% check input
if isempty(spiketimes) || isempty(cids)
    error('Input arguments missing. Extract spikes first.')
elseif isempty(Srise)
    error('Input arguments missing. Define Srise first.')
elseif isempty(stimuli_parameters)
    error('Stimuli parameters file missing.')
end

% check if all trials have TTL
NStim =  size(stimuli_parameters.Stm, 1);
if NStim ~= size(Srise)
    error('Length Srise and NStim do not correspond')
else
    disp(['Stimuli in session: ' num2str(NStim)])
end

% select correct analysis window
if strcmp(stimuli_parameters.Par.Rec, 'SOM')
    PreT = str2double(stimuli_parameters.Par.SomatosensoryISI)/4; % amount of msec. to include before Srise;
    PostT = (str2double(stimuli_parameters.Par.SomatosensoryStimTime) + str2double(stimuli_parameters.Par.SomatosensoryISI)/4); % amount of msec. to include after Sfall;
    % -- NEW -- %
elseif strcmp(stimuli_parameters.Par.Rec, 'SxA') % TO TEST
    if length(unique(~isnan(stimuli_parameters.Stm.SomAudSOA))) > 2 % multiple SOA delays
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
elseif strcmp(stimuli_parameters.Par.Rec, 'FRA') || strcmp(stimuli_parameters.Par.Rec, 'Opt') || strcmp(stimuli_parameters.Par.Rec, 'OptoFRA')
    PreT  = str2double(stimuli_parameters.Par.FRAStimTime);
    PostT = str2double(stimuli_parameters.Par.FRAPostTime);
elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
    PreT = (str2double(stimuli_parameters.Par.AMPreTime) + (str2double(stimuli_parameters.Par.AMPostTime)/4));
    PostT = (str2double(stimuli_parameters.Par.AMStimTime) + (str2double(stimuli_parameters.Par.AMPostTime)/2));
end


% initialize variables (does nothing currently)
SpkT = [];
%aligned_spikes = {}; % nStim x Ncids

% actual spike alignment done here
for cluster = 1:cids

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

%aligned_spikes = [aligned_spikes; SpkT]; % stimuli x units recording
end
