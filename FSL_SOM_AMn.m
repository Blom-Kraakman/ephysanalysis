function [SOM] = FSL_SOM_AMn(stimuli_parameters, aligned_spikes, cids)
% give median first spike latenct following noise or somatosensory
% stimulation

% stimulus parameters
if strcmp(stimuli_parameters.Par.Rec, 'SOM')
    x = stimuli_parameters.Stm.Amplitude;
elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
    x = [stimuli_parameters.Stm.Intensity];
    x((x ~= 15) & (x ~= 30) & (x ~= 45) & (x ~= 60)) = 0;
end

UAmp = unique(x);
NAmp = length(UAmp);
NClu = length(cids); % cluster info

% analysis period
%winStart = 0;
%winEnd   = 0.3; % changed from FRA

% initiation variables
NTrials = nan(NAmp, NClu); % number of trials
MedFSL = nan(NAmp, NClu); % median first spike latency
IqrFSL = nan(NAmp, NClu); % IQR first spike latency

for cluster = 1:NClu

    % Spike count analysis
    for amplitude = 1:NAmp
            sel = x == UAmp(amplitude);
            NTrials(amplitude,cluster) = sum(sel);

            if (NTrials(amplitude,cluster) == 0); continue; end

            %count spikes
            tempSpiketimes = aligned_spikes(sel,cluster);

            % SCnt = nan(NTrials(amplitude,cluster), 1);
            % for t = 1:length(tempSpiketimes)
            %     if (isnan(tempSpiketimes{t})); continue; end
            %     if(isempty(tempSpiketimes{t}))
            %         SCnt(t) = 0;
            %     else
            %         SCnt(t) = sum(tempSpiketimes{t} > winStart & tempSpiketimes{t} < winEnd);
            %     end
            % end

            %FSL
            fsl = inf(length(tempSpiketimes), 1);
            for t = 1:length(tempSpiketimes)
                if (isnan(tempSpiketimes{t}))
                    continue
                end

                spks = tempSpiketimes{t};
                spks = spks (spks > 0);

                if (~isempty(spks))
                    fsl(t) = min(spks);
                end
            end
            
            MedFSL(amplitude,cluster) = median(fsl);
            IqrFSL(amplitude,cluster) = iqr(fsl);

    end % amplitude loop
    clearvars('tempSpiketimes');

end % cluster loop


% Output data Experiment meta data

% cluster
SOM.NClu = NClu;

% stimulus parameters
SOM.UAmp = UAmp;
SOM.NAmp = NAmp;

% results - number of trials
SOM.NTrials = NTrials;

% results - spike latency
SOM.MedFSL = MedFSL;
SOM.IqrFSL = IqrFSL;

%% Plot FSL data
% plot medFSL
f = plot(SOM.UAmp, SOM.MedFSL, '-x');

% format figure
legend()
ylabel('FSL (s)');
xticks(UAmp);
xlim([(UAmp(1) - 2), (UAmp(end) +2)]);

if strcmp(stimuli_parameters.Par.Rec, 'SOM')
    xlabel('Condition');
    xticklabels({'stim off' 'stim on'})
elseif strcmp(stimuli_parameters.Par.Rec, 'AMn')
    xlabel('Codition (dB SPL)');
end

end