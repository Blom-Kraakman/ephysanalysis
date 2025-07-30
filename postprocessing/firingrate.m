function [FR, nrSpikes] = firingrate(SpkT, PreT, PostT)
% calculate firing rate in defined time window
% INPUT - SpkT: cell array with spike times aligned to stimulus onset
% PreT & PostT (in sec) arrays

% input check
if ~isequal(size(SpkT,1), length(PreT), length(PostT))
    error('Input does not contain same amount of trials and analysis windows')
end

% initiate variables
nClust = size(SpkT,2);
nTrials = size(SpkT,1);
FR = nan(nTrials, nClust);
nrSpikes = nan(nTrials, nClust);
window = PostT - PreT;

% firing rate calculation
for cluster = 1:nClust
    for trial = 1:nTrials
        
        % get spike times for 1 trial
        S = SpkT{trial, cluster};
        
        % count spikes in time windows
        tnrSpikes = sum(S >= -PreT(trial) & S <= PostT(trial));

        % convert to firing rate (Hz)
        nrSpikes(trial, cluster) = tnrSpikes;
        FR(trial, cluster) = tnrSpikes/abs(window(trial));

    end
end
end