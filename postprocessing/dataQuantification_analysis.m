function dataQuantification_analysis(animal, OutPath, BehaviorPath, session)
% Load relevant data files
[cids, stimuli_parameters, aligned_spikes, ~, ~, ~, onsetDelay, StimResponseFiring] = loadData(OutPath, session, BehaviorPath);

%% 1. calculate firing rates

% Initiate table
MouseNum = repmat(animal,[length(cids'),1]);
%Session = repmat(session,[length(cids'),1]);
unitResponses = table(MouseNum, cids');
unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster'};

% define analysis window - entire stimulus period
if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') % && strcmp(stimuli_parameters.Par.Rec, "SxA")
    PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
    PostT = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
    %PostT = 0.25; % captures initial noise period & half of vibrotac (in noise) period
    %PostT = 0.5; % whole vibrotac + dual mode period
elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square') % && strcmp(stimuli_parameters.Par.Rec, "SOM")
    if strcmp(stimuli_parameters.Par.Rec, 'SxA') && length(str2num(stimuli_parameters.Par.SomAudSOA)) > 2
        PreT = 0.1;
        PostT = 0.13;
    else
        PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
        PostT = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
        %PostT = 0.1; % best way to capture onset stimulus
    end
elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'BiSine')
    PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
    PostT = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
    if max(stimuli_parameters.Stm.Md)
        PreT = str2double(stimuli_parameters.Par.AMPreTime)/1000;
        PostT = str2double(stimuli_parameters.Par.AMPostTime)/1000;
    else
        %PreT = (str2double(stimuli_parameters.Par.AMISI)/3)/1000; % baseline period
        PreT = 0.2;
        PostT = 0.2; % limited by previous recordings
    end
elseif strcmp(stimuli_parameters.Par.Rec, "FRA")
    PreT = 0.2;
    PostT = (str2double(stimuli_parameters.Par.FRAStimTime)/1000);
else
    error("No analysis window defined")
end

% calculate firing rate (Hz) in time window
%[baselineRate, stimulusRate] = firingrate(aligned_spikes, PreT, PostT, onsetDelay);

% format
PostT = repmat(PostT, size(aligned_spikes,1),1);
PreT = repmat(PreT, size(aligned_spikes,1),1);

% calculate baseline FR
baselineRate = firingrate(aligned_spikes, PreT, zeros(size(aligned_spikes,1), 1));

% calculate stimulus induced FR
stimulusRate = firingrate(aligned_spikes, onsetDelay, (PostT+onsetDelay));

% calculate stimulus induced change in firing rate
dstimulusRate = stimulusRate - baselineRate;

clear MouseNum

%% 2. calculate mean stimulus induced change in firing & first spike latency

% set parameter space
%[uAmp, uFreq, conditions] = selectparameters(stimuli_parameters);

% select correct amplitude and frequency parameters
if strcmp(stimuli_parameters.Par.Rec, "AMn") % noise session
    uparamA = unique(stimuli_parameters.Stm.Intensity);
    uparamB = unique(stimuli_parameters.Stm.Mf);
elseif strcmp(stimuli_parameters.Par.Rec, "FRA")
    uparamA = unique(stimuli_parameters.Stm.Intensity);
    uparamB = unique(stimuli_parameters.Stm.Freq);
elseif strcmp(stimuli_parameters.Par.Rec, "SxA") && strcmp(stimuli_parameters.Par.SomatosensoryWaveform, "Square") % SxA pressure session
    uparamA = unique(stimuli_parameters.Stm.Amplitude);
    uparamB = unique(stimuli_parameters.Stm.AudIntensity); % sound variable
else % SxA vibrotac and SOM sessions
    uparamA = unique(stimuli_parameters.Stm.Amplitude);
    uparamB =  unique(stimuli_parameters.Stm.SomFreq);
end

% select multiple stimuli types in one session
if strcmp(stimuli_parameters.Par.Rec, "SxA") %dual modal session
    conditions = {'OO', 'OA', 'SO', 'SA'};
else % unimodal sessions
    conditions = 1;
end

%nAmp = length(uAmp);
%nFreq = length(uFreq);
nparamA = length(uparamA);
nparamB = length(uparamB);
NClu = length(cids);

% initialize output variables
firing_mean = nan(nparamB, nparamA, length(conditions), NClu);
MedFSL = nan(nparamB, nparamA, length(conditions), NClu); %nan(nAmp, nFreq, NClu); % median first spike latency
IqrFSL = nan(nparamB, nparamA, length(conditions), NClu); %nan(nAmp, nFreq, NClu); % IQR first spike latency

% mean dfiring rate per stimulus combination for all units
for condition = 1:length(conditions)
    for freq = 1:nparamB
        for amp = 1:nparamA
            if strcmp(stimuli_parameters.Par.Rec, "SxA")
                if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
                    index = stimuli_parameters.Stm.Amplitude == uparamA(amp) & stimuli_parameters.Stm.AudIntensity == uparamB(freq) & strcmp(stimuli_parameters.Stm.MMType, conditions{condition});
                elseif strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
                    index = stimuli_parameters.Stm.Amplitude == uparamA(amp) & stimuli_parameters.Stm.SomFreq == uparamB(freq) & strcmp(stimuli_parameters.Stm.MMType, conditions{condition});
                end
            elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
                index = stimuli_parameters.Stm.Intensity == uparamA(amp) & stimuli_parameters.Stm.Mf == uparamB(freq);
            elseif strcmp(stimuli_parameters.Par.Rec, "FRA")
                index = stimuli_parameters.Stm.Intensity == uparamA(amp) & stimuli_parameters.Stm.Freq == uparamB(freq);
            elseif strcmp(stimuli_parameters.Par.Rec, "SOM")
                index = stimuli_parameters.Stm.Amplitude == uparamA(amp) & stimuli_parameters.Stm.SomFreq == uparamB(freq);
            end

            if sum(index) == 0; continue; end

            % firing rate analysis
            firing_mean(freq, amp, condition, :) = mean(dstimulusRate(index, 1:NClu));

            % first spike time analysis
            for cluster = 1:size(aligned_spikes,2)
                tempSpiketimes = aligned_spikes(index,cluster);
                fsl = inf(size(tempSpiketimes,1),1);

                for trial = 1:size(tempSpiketimes,1)
                    if (isnan(tempSpiketimes{trial}))
                        continue
                    end

                    spks = tempSpiketimes{trial};
                    spks = spks (spks > 0);

                    if (~isempty(spks))
                        fsl(trial) = min(spks);
                    end
                end

                MedFSL(freq,amp,condition,cluster) = median(fsl);
                IqrFSL(freq,amp,condition,cluster) = iqr(fsl);

            end
        end
    end
end

clear freq amp condition cluster trial tempSpiketimes index spks

%% 3. quantify reactive units
% ! edit AMn indexing if needed (needed in earlier data sets)

hval = nan(nparamB,nparamA,NClu);
pval = nan(nparamB,nparamA,NClu);

for cond = 1:length(conditions)
   
    [responsive, thval, tpval] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, conditions(cond), uparamA, uparamB);

    % fill hval and pval matrix for all stim combinations
    if strcmp(conditions(cond), 'OO')
        hval(1,1,:) = thval;
        pval(1,1,:) = tpval;
    elseif strcmp(conditions(cond), 'OA')
        hval(2:end,1,:) = thval(2:end,1,:);
        pval(2:end,1,:) = tpval(2:end,1,:);
    elseif strcmp(conditions(cond), 'SO')
        hval(1,2:end,:) = thval(1,2:end,:);
        pval(1,2:end,:) = tpval(1,2:end,:);
    elseif strcmp(conditions(cond), 'SA')
        hval(2:end,2:end,:) = thval(2:end,2:end,:);
        pval(2:end,2:end,:) = tpval(2:end,2:end,:);
    end
            % index responsive units
    idx = max(responsive == unitResponses.Cluster, [], 2);

    if isempty(idx)
        idx = false(NClu,1);
    end

    % add responsive true/false per cluster per condition
    unitResponses = addvars(unitResponses, idx);
end

% add OA sound steady state
if strcmp(stimuli_parameters.Par.Rec, "SxA") && str2num(stimuli_parameters.Par.SomAudSOA) > 0

    onsetDelay_OA = repmat(max(onsetDelay), length(onsetDelay), 1);

    % calculate stimulus induced FR, time window now same as som
    stimulusRate = firingrate(aligned_spikes, onsetDelay_OA, (PostT+onsetDelay_OA));

    responsive = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, 'OA', uparamA, uparamB); % units responsive

    % index responsive units
    idx = max(responsive == unitResponses.Cluster, [], 2);

    if isempty(idx)
        idx = false(NClu,1);
    end

    % add responsive true/false per cluster per condition
    unitResponses = addvars(unitResponses, idx);

    % modify table comlumn names and add to struct
    unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster', 'OO', 'OA', 'SO', 'SA', 'OA+0.25'};

elseif strcmp(stimuli_parameters.Par.Rec, "SOM")
    unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster', 'SOM'};

elseif strcmp(stimuli_parameters.Par.Rec, "SxA") && str2num(stimuli_parameters.Par.SomAudSOA) == 0
    unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster', 'OO', 'OA', 'SO', 'SA'};

elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
    unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster', 'AM'};
end

StimResponseFiring.unitResponses = unitResponses;
StimResponseFiring.MouseNum = stimuli_parameters.Par.MouseNum;
StimResponseFiring.session = stimuli_parameters.Par.Set;
StimResponseFiring.type = stimuli_parameters.Par.Rec;
StimResponseFiring.conditions = conditions;
StimResponseFiring.cids = cids;
StimResponseFiring.amplitudes = uparamA;
StimResponseFiring.frequencies = uparamB;
StimResponseFiring.PreT = PreT;
StimResponseFiring.PostT = PostT;
StimResponseFiring.onsetDelay = onsetDelay;
StimResponseFiring.baselineRate = baselineRate;
StimResponseFiring.stimulusRate = stimulusRate;
StimResponseFiring.dfiringRate = dstimulusRate;
StimResponseFiring.firing_mean = firing_mean;
StimResponseFiring.FSLmed = MedFSL;
StimResponseFiring.FSLiqr = IqrFSL;
StimResponseFiring.pvalue = pval;
StimResponseFiring.hvalue = hval;

clear cond index idx

%% 3.5 Save
filename = sprintf('M%.2i_S%.2i_%s_ResponseProperties', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
save(fullfile([OutPath '\ResponseProperties'], filename), "StimResponseFiring")
    
end
