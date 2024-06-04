function [results, resp_cids] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids)

% quantify responsive units
% sig diff firing rate between stim period vs baseline
% INPUT: session parameters, cids, firing rate (trials x units) during specific session
% OUTPUT: struct containing id of responsive units

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

% initiate struct
%results.session = stimuli_parameters.Par.Set;
%results.type = stimuli_parameters.Par.Rec;
%results.cids = cids;
%results.ampitudes = uAmp;
%results.frequencies = uFreq;
results.pvalue = nan(nAmp, nFreq, nClusters);
results.zval = nan(nAmp, nFreq, nClusters);
results.ranksum = nan(nAmp, nFreq, nClusters);
results.hvalue = nan(nAmp, nFreq, nClusters);
results.responsive = [];

% stats test difference baseline vs stimulus period
for cluster = 1:nClusters

    %signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    for freq = 1:nFreq
        for condition = 1:nAmp

            % session specific parameters
            if strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) == 0 % noise
                index = tAmp == uAmp(condition);
                control = tAmp == uAmp(1);
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            elseif strcmp(stimuli_parameters.Par.Rec, 'AMn') && max(stimuli_parameters.Stm.Mf) ~= 0 % AM
                index = stimuli_parameters.Stm.Mf == uFreq(condition);
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
                index = (strcmp(stimuli_parameters.Stm.MMType, 'OA')) & (stimuli_parameters.Stm.Amplitude == uAmp(condition)) & (stimuli_parameters.Stm.SomFreq == uFreq(freq));
                baseline = stimulusRate(control,cluster);
                stim_evoked = stimulusRate(index,cluster);
            end

            if isempty(stim_evoked) || isempty(baseline)
                continue
            end

            % stat test
            [p,h,stats] = ranksum(baseline, stim_evoked, 'alpha', 0.01); % nonpaired version: experimental diff from control
            %[p,h,stats] = signrank(baseline, stim_evoked, 'alpha', 0.01); % exp trial diff from baseline
            results.pvalue(condition, freq, cluster) = p;
            results.zval(condition, freq, cluster) = stats.zval;
            results.ranksum(condition, freq, cluster) = stats.ranksum;
            results.hvalue(condition, freq, cluster) = h;
        end
    end

    % % unit responsive if sig diff for at least one condition
    % if (strcmp(stimuli_parameters.Par.Rec, 'SxA')) && (max(max(results.hvalue(2:nAmp, 2:nFreq, cluster))) == 1) % && (results.hvalue(1, 1, cluster) == 0)
    %     results.responsive = [results.responsive, results.cids(cluster)];
    % elseif (strcmp(stimuli_parameters.Par.Rec, 'AMn')) && max(results.hvalue(:, :, cluster)) % && (results.hvalue(1, 1, cluster) == 0))
    %     results.responsive = [results.responsive, results.cids(cluster)];
    % elseif (strcmp(stimuli_parameters.Par.Rec, 'SOM'))
    %     if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') && max(results.hvalue(:, :, cluster))
    %         results.responsive = [results.responsive, results.cids(cluster)];
    %     elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') && max(max(results.hvalue(:, :, cluster)))
    %         results.responsive = [results.responsive, results.cids(cluster)];
    %     end
    % end
end

% unit responsive if sig diff for at least one condition
if strcmp(stimuli_parameters.Par.Rec, 'AMn')
    index = squeeze(max(results.hvalue)) == 1;
elseif (strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')) && (strcmp(stimuli_parameters.Par.Rec, 'SOM'))
    index = squeeze(max(results.hvalue)) == 1;
else
    index = squeeze(max(max(results.hvalue))) == 1;
end
results.responsive = results.cids(index);

%responsive units contains all
resp_cids = find(ismember(results.cids, results.responsive)); % position of responsive units in cids

end