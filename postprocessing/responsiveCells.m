function [responsive, hval, pval] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, condition, uparamA, uparamB)
% quantify responsive units
% sig diff firing rate between stim period vs baseline
% INPUT: session parameters, cids, firing rate (trials x units) during specific session
% OUTPUT: struct containing id of responsive units

if ~exist('condition', 'var')
    condition = [];
end

if ~exist('uparamA', 'var') && strcmp(condition, 'OA')
    uparamA = [];
end

nparamA = length(uparamA);
nparamB = length(uparamB);
NClu = length(cids);

% intiate variables
pval = nan(nparamB, nparamA, NClu);
zval = nan(nparamB, nparamA, NClu);
hval = nan(nparamB, nparamA, NClu);

% set alpha value
alpha_val =  0.05/(nparamB*nparamA);

% stats test difference baseline vs stimulus period

%signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
for freq = 1:nparamB
    for amp = 1:nparamA

        % select session specific parameters
        [baseline, stim_evoked] = getTrials(stimuli_parameters, baselineRate, stimulusRate, amp, freq, uparamA, uparamB, condition);

        if isempty(stim_evoked) || isempty(baseline)
            continue
        end

        % stat test
        for cluster = 1:NClu
            [p,h,stats] = ranksum(baseline(:, cluster), stim_evoked(:, cluster), 'alpha',  alpha_val); % nonpaired version: experimental diff from control
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



%% LOCAL FUNCTIONS
% get the firing rate per stim combination
    function [baseline, stim_evoked] = getTrials(stimuli_parameters, baselineRate, stimulusRate, amp, freq, uparamA, uparamB, condition)
        baseline = [];
        stim_evoked = [];

        switch stimuli_parameters.Par.Rec
            case 'AMn'
                if max(stimuli_parameters.Stm.Mf) == 0 % bbnoise
                    index = stimuli_parameters.Stm.Intensity == uparamA(amp);
                    control = stimuli_parameters.Stm.Intensity == uparamA(1);
                    baseline = stimulusRate(control,:);
                    stim_evoked = stimulusRate(index,:);
                else % AM noise
                    index = stimuli_parameters.Stm.Mf == uparamB(freq);
                    baseline = baselineRate(index,:);
                    stim_evoked = stimulusRate(index,:);
                end

            case 'FRA'
                index = (stimuli_parameters.Stm.Intensity == uparamA(amp)) & (stimuli_parameters.Stm.Freq == uparamB(freq));
                control = stimuli_parameters.Stm.Freq == 0;
                baseline = stimulusRate(control,:);
                stim_evoked = stimulusRate(index,:);

            case 'SOM'
                waveform = stimuli_parameters.Par.SomatosensoryWaveform;
                switch waveform
                    case 'UniSine' % vibrotac
                        index = (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.SomFreq == uparamB(freq));
                        control = (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.SomFreq == 0);
                    case 'Square' % pressure
                        index = stimuli_parameters.Stm.Amplitude == uparamA(amp);
                        control = stimuli_parameters.Stm.Amplitude == 0;
                    case 'BiSine' % vibrotac (piezo)
                        index = (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.SomFreq == uparamB(freq));
                        control = (stimuli_parameters.Stm.Amplitude == 0) & (stimuli_parameters.Stm.SomFreq == 0);
                end
                baseline = stimulusRate(control,:);
                stim_evoked = stimulusRate(index,:);

            case 'SxA' % multimodal
                control = strcmp(stimuli_parameters.Stm.MMType, 'OO');
                waveform = stimuli_parameters.Par.SomatosensoryWaveform;
                if strcmp(waveform, 'Square') && strcmp(condition, 'OA') % sound only
                    index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.AudIntensity == uparamB(freq));
                elseif strcmp(waveform, 'Square') && (strcmp(condition, 'SA') || strcmp(condition, 'SO')) % multimodal or som only
                    index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.AudIntensity == uparamB(freq));
                elseif strcmp(waveform, 'UniSine')
                    index = (strcmp(stimuli_parameters.Stm.MMType, condition)) & (stimuli_parameters.Stm.Amplitude == uparamA(amp)) & (stimuli_parameters.Stm.SomFreq == uparamB(freq));
                end
                baseline = stimulusRate(control,:);

        end
    end
end
