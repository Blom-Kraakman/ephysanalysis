function [responsive_units, resp_cids] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, session)

% quantify responsive units
% sig diff firing rate between stim period vs baseline
% INPUT: session parameters, cids, firing rate (trials x units) during specific session
% OUTPUT: struct containing id of responsive units

% set parameters
if strcmp(stimuli_parameters.Par.Rec, 'AMn')
    uFreq = unique(stimuli_parameters.Stm.Mf);
    uAmp = unique(stimuli_parameters.Stm.Intensity); %intensity for M4
else
    uAmp = unique(stimuli_parameters.Stm.Amplitude);
    uFreq = unique(stimuli_parameters.Stm.SomFreq);
end

nAmp = length(uAmp);
nFreq = length(uFreq);

nClusters = length(cids);

% initiate struct
results.session = stimuli_parameters.Par.Set;
results.type = stimuli_parameters.Par.Rec;
results.cids = cids;
results.ampitudes = uAmp;
results.frequencies = uFreq;
results.pvalue = nan(nAmp, nFreq, nClusters);
results.signrank = nan(nAmp, nFreq, nClusters);
results.hvalue = nan(nAmp, nFreq, nClusters);
results.responsive = [];

% stats test difference baseline vs stimulus period
for cluster = 1:nClusters

    %signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    for freq = 1:nFreq
        for condition = 1:nAmp

            % session specific parameters
            if strcmp(stimuli_parameters.Par.Rec, 'AMn') % noise & AM
                %index = stimuli_parameters.Stm.Mf == uFreq(condition);
                %m4 index
                index = stimuli_parameters.Stm.Intensity == uAmp(condition);
                baseline = baselineRate(index,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'SOM')
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') %vibrotac
                control = (stimuli_parameters.Stm.Amplitude == uAmp(condition)) & (stimuli_parameters.Stm.SomFreq == 0);
                index = (stimuli_parameters.Stm.Amplitude == uAmp(condition)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') % pressure
                    control = stimuli_parameters.Stm.Amplitude == 0;
                    index = stimuli_parameters.Stm.Amplitude == uAmp(condition);
                    baseline = stimulusRate(control,cluster);
                    stim_evoked = stimulusRate(index,cluster);
                end
            elseif strcmp(stimuli_parameters.Par.Rec, 'SxA') % multimodal
                control = strcmp(stimuli_parameters.Stm.MMType, 'OO');
                %control = (stimuli_parameters.Stm.Amplitude == 0) & (stimuli_parameters.Stm.SomFreq == 0);
                index = (strcmp(stimuli_parameters.Stm.MMType, 'SO')) & (stimuli_parameters.Stm.Amplitude == uAmp(condition)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            end

            %[p,h,stats] = signrank(stimulusRate(control,cluster), stimulusRate(index,cluster), 'alpha', 0.01); % different from control trials?
            % stat test
            [p,h,stats] = signrank(baseline, stim_evoked, 'alpha', 0.01); % different from control trials?
            results.pvalue(condition, freq, cluster) = p;
            results.signrank(condition, freq, cluster) = stats.signedrank;
            results.hvalue(condition, freq, cluster) = h;
        end
    end

    % unit responsive if sig diff for at least one condition
    if (strcmp(stimuli_parameters.Par.Rec, 'SxA')) && (max(max(results.hvalue(2:nAmp, 2:nFreq, cluster))) == 1) % && (results.hvalue(1, 1, cluster) == 0)
        results.responsive = [results.responsive, results.cids(cluster)];
    elseif (strcmp(stimuli_parameters.Par.Rec, 'AMn')) && max(results.hvalue(:, :, cluster)) % && (results.hvalue(1, 1, cluster) == 0))
        results.responsive = [results.responsive, results.cids(cluster)];
    elseif (strcmp(stimuli_parameters.Par.Rec, 'SOM'))
        if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') && max(results.hvalue(:, :, cluster))
            results.responsive = [results.responsive, results.cids(cluster)];
        elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') && max(max(results.hvalue(:, :, cluster)))
            results.responsive = [results.responsive, results.cids(cluster)];
        end
    end
end

%responsive units contains all
responsive_units(str2double(stimuli_parameters.Par.Set)) = results;
resp_cids = find(ismember(responsive_units(session).cids, responsive_units(session).responsive)); % position of responsive units in cids

end