function [baseFR, expFR, baseSpikes, expSpikes] = firingrate(SpkT, PreT, PostT, StimOn)
% calculate baseline firing rate
% SpkT: cell array with spike times aligned to stimulus onset
nClust = size(SpkT,2);
nTrials = size(SpkT,1);
baseFR = nan(nTrials, nClust);
expFR = nan(nTrials, nClust);

for cluster = 1:size(SpkT,2)
    for trial=1:size(SpkT,1)
        
        % get spike times per trial
        S = SpkT{trial, cluster};
        
        % count spikes in time windows
        baseSpikes = sum(S >= -PreT & S <= 0); % preT - 0
        expSpikes = sum(S > StimOn(trial) & S <= (StimOn(trial) + PostT)); % StimOn moment - postT

        % convert to firing rate (Hz)
        baseFR(trial, cluster) = baseSpikes/abs(PreT);
        expFR(trial, cluster) = expSpikes/abs(PostT);
    end
end
end