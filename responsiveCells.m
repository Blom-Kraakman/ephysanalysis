function [responsive, hval, pval] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, condition, nparamA, uparamA,nparamB, uparamB)

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

% % set parameters
% [uparamA, uparamB, ~] = selectparameters(stimuli_parameters);
% nparamA = length(uparamA);
% nparamB = length(uparamB); % contains aud amp in SxA pressure sessions
NClu = length(cids);

pval = nan(nparamB, nparamA, NClu);
zval = nan(nparamB, nparamA, NClu);
hval = nan(nparamB, nparamA, NClu);

% stats test difference baseline vs stimulus period
for cluster = 1:NClu

    %signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    for freq = 1:nparamB
        for amp = 1:nparamA
            % session specific parameters
            if strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) == 0 % bbnoise
                index = tAmp == uparamA(amp);
                control = tAmp == uparamA(1);
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) ~= 0 % AM noise
                index = stimuli_parameters.Stm.Mf == uparamB(amp);
                baseline = baselineRate(index,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'FRA')
                index = (stimuli_parameters.Stm.Intensity == uparamA(amp)) & (stimuli_parameters.Stm.Freq == uparamB(freq));
                control = stimuli_parameters.Stm.Freq == 0;
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') %vibrotac
                    index = (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.SomFreq == uparamB(freq));
                    control = (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.SomFreq == 0);
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') % pressure
                    index = stimuli_parameters.Stm.Amplitude == uparamA(amp);
                    control = stimuli_parameters.Stm.Amplitude == 0;
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                end
            elseif strcmp(stimuli_parameters.Par.Rec, 'SxA') % multimodal, vibrotac
                control = strcmp(stimuli_parameters.Stm.MMType, 'OO');
                %control = (stimuli_parameters.Stm.Amplitude == 0) & (stimuli_parameters.Stm.SomFreq == 0);
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
                    index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.AudIntensity == uparamB(freq));
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
                    index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.SomFreq == uparamB(freq));
                end
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            end

            if isempty(stim_evoked) || isempty(baseline)
                continue
            end

            % stat test
            alpha_val =  0.05/(nparamB*nparamA);
            [p,h,stats] = ranksum(baseline, stim_evoked, 'alpha',  alpha_val); % nonpaired version: experimental diff from control
            %[p,h,stats] = signrank(baseline, stim_evoked, 'alpha', 0.01); % paired versiom: exp trial diff from baseline
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
