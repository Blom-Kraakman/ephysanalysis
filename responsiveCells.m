function [responsive] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, condition)

% quantify responsive units
% sig diff firing rate between stim period vs baseline
% INPUT: session parameters, cids, firing rate (trials x units) during specific session
% OUTPUT: struct containing id of responsive units

if ~exist('condition', 'var')
    condition = [];
end

% set parameters
if strcmp(stimuli_parameters.Par.Rec, 'AMn') % noise AM
    tAmp = stimuli_parameters.Stm.Intensity;
    uAmp = unique(tAmp);
    uFreq = unique(stimuli_parameters.Stm.Mf);
else
    tAmp = stimuli_parameters.Stm.Amplitude;
    uAmp = unique(tAmp);
    uFreq = unique(stimuli_parameters.Stm.SomFreq);
end

nAmp = length(uAmp);
nFreq = length(uFreq);
nClusters = length(cids);

pval = nan(nAmp, nFreq, nClusters);
zval = nan(nAmp, nFreq, nClusters);
hval = nan(nAmp, nFreq, nClusters);

% stats test difference baseline vs stimulus period
for cluster = 1:nClusters

    %signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    for freq = 1:nFreq
        for amp = 1:nAmp

            % session specific parameters
            if strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) == 0 % noise
                index = tAmp == uAmp(amp);
                control = tAmp == uAmp(1);
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) ~= 0 % AM
                index = stimuli_parameters.Stm.Mf == uFreq(amp);
                baseline = baselineRate(index,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') %vibrotac
                    control = (stimuli_parameters.Stm.Amplitude == uAmp(amp)) & (stimuli_parameters.Stm.SomFreq == 0);
                    index = (stimuli_parameters.Stm.Amplitude == uAmp(amp)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') % pressure
                    control = stimuli_parameters.Stm.Amplitude == 0;
                    index = stimuli_parameters.Stm.Amplitude == uAmp(amp);
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                end
            elseif strcmp(stimuli_parameters.Par.Rec, 'SxA') % multimodal
                control = strcmp(stimuli_parameters.Stm.MMType, 'OO');
                %control = (stimuli_parameters.Stm.Amplitude == 0) & (stimuli_parameters.Stm.SomFreq == 0);
                index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.Amplitude == uAmp(amp)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            end

            if isempty(stim_evoked) || isempty(baseline)
                continue
            end

            % stat test
            [p,h,stats] = ranksum(baseline, stim_evoked, 'alpha', 0.01); % nonpaired version: experimental diff from control
            %[p,h,stats] = signrank(baseline, stim_evoked, 'alpha', 0.01); % exp trial diff from baseline
            pval(amp, freq, cluster) = p;
            zval(amp, freq, cluster) = stats.zval;
            hval(amp, freq, cluster) = h;
        end
    end

end

% unit responsive if sig diff for at least one condition
if strcmp(stimuli_parameters.Par.Rec, 'AMn')
    index = squeeze(max(hval)) == 1;
elseif (strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')) && (strcmp(stimuli_parameters.Par.Rec, 'SOM'))
    index = squeeze(max(hval)) == 1;
else
    index = squeeze(max(max(hval))) == 1;
end

% responsive units
responsive = cids(index);

end
