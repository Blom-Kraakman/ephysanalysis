function [aligned_spikes] = alignspikes(spiketimes, cids, Srise, NStim, PreT, PostT, Fs)
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
elseif isempty(PreT) || isempty(PostT)
    error('Input arguments missing. Define alignement window.')
end

% check if all trials have TTL
if NStim ~= size(Srise)
    error('Length Srise and NStim do not correspond')
end

% initialize variables (does nothing currently)
SpkT = [];
aligned_spikes = {}; % nStim x Ncids

% actual spike alignment done here
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

aligned_spikes = [aligned_spikes; SpkT]; % stimuli x units recording
fprintf('spike alignment done\n');

end
