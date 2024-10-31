function [responsive, hval, pval] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, condition)

% quantify responsive units
% sig diff firing rate between stim period vs baseline
% INPUT: session parameters, cids, firing rate (trials x units) during specific session
% OUTPUT: struct containing id of responsive units

if ~exist('condition', 'var')
    condition = [];
end

% % set parameters
% if strcmp(stimuli_parameters.Par.Rec, 'AMn') % noise AM
%     tAmp = stimuli_parameters.Stm.Intensity;
%     uFreq = unique(stimuli_parameters.Stm.Mf);
% elseif strcmp(stimuli_parameters.Par.Rec, 'FRA')
%     tAmp = stimuli_parameters.Stm.Intensity;
%     uFreq = unique(stimuli_parameters.Stm.Freq);
% else
%     tAmp = stimuli_parameters.Stm.Amplitude;
%     uFreq = unique(stimuli_parameters.Stm.SomFreq);
% end
% 
% uAmp = unique(tAmp);
% nAmp = length(uAmp);
% nFreq = length(uFreq);
% NClu = length(cids);

% set parameters
[uAmp, uFreq, ~] = selectparameters(stimuli_parameters);
nAmp = length(uAmp);
nFreq = length(uFreq); % contains aud amp in SxA pressure sessions
NClu = length(cids);

pval = nan(nFreq, nAmp, NClu);
zval = nan(nFreq, nAmp, NClu);
hval = nan(nFreq, nAmp, NClu);

% stats test difference baseline vs stimulus period
for cluster = 1:NClu

    %signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    for freq = 1:nFreq
        for amp = 1:nAmp
            % session specific parameters
            if strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) == 0 % bbnoise
                index = tAmp == uAmp(amp);
                control = tAmp == uAmp(1);
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) ~= 0 % AM noise
                index = stimuli_parameters.Stm.Mf == uFreq(amp);
                baseline = baselineRate(index,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'FRA')
                index = (stimuli_parameters.Stm.Intensity == uAmp(amp)) & (stimuli_parameters.Stm.Freq == uFreq(freq));
                control = stimuli_parameters.Stm.Freq == 0;
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') %vibrotac
                    index = (stimuli_parameters.Stm.Amplitude == uAmp(amp)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                    control = (stimuli_parameters.Stm.Amplitude == uAmp(amp)) & (stimuli_parameters.Stm.SomFreq == 0);
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') % pressure
                    index = stimuli_parameters.Stm.Amplitude == uAmp(amp);
                    control = stimuli_parameters.Stm.Amplitude == 0;
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                end
            elseif strcmp(stimuli_parameters.Par.Rec, 'SxA') % multimodal, vibrotac
                control = strcmp(stimuli_parameters.Stm.MMType, 'OO');
                %control = (stimuli_parameters.Stm.Amplitude == 0) & (stimuli_parameters.Stm.SomFreq == 0);
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
                    index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.Amplitude == uAmp(amp)) & (stimuli_parameters.Stm.AudIntensity == uFreq(freq));
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
                    index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.Amplitude == uAmp(amp)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                end
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            end

            if isempty(stim_evoked) || isempty(baseline)
                continue
            end

            % stat test
            %[p,h,stats] = ranksum(baseline, stim_evoked, 'alpha', 0.05/(nAmp*nFreq)); % stricter alpha correction
            [p,h,stats] = ranksum(baseline, stim_evoked, 'alpha', 0.01); % nonpaired version: experimental diff from control
            %[p,h,stats] = signrank(baseline, stim_evoked, 'alpha', 0.01); % exp trial diff from baseline
            pval(freq, amp,cluster) = p;
            zval(freq, amp, cluster) = stats.zval;
            hval(freq, amp,cluster) = h;
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
