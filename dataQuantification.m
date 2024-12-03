%% Quantify results
% To be done after PostCuration_BK.m
% General data input: single units (array) stimuli order (table), aligned spikes (cell array), 
% Post processing data steps:
%   1. get firing rate during set time frame for each trial of each unit
%   2. quatify which cluster responds to which event
%       responsive unit = sig different firing after stim vs baseline
%   3. quantify response strength, make various plots
%   4. combine data from multiple animals

clearvars

%SxA, SOM
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
%session = [4, 6]; % M9

% Data paths to sorted data, behavioral data file, output folder

animal = 16;
session = 6;

BehaviorPath = ['D:\DATA\Behavioral Stimuli\M' num2str(animal, 2)]; % stimuli parameters
OutPath = ['D:\DATA\Processed\M' num2str(animal) '\ICX\rec2']; % output directory

[cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, onsetDelay, StimResponseFiring] = loadData(OutPath, session, BehaviorPath);

%% 1. calculate firing rates

% Initiate table
MouseNum = repmat(animal,[length(cids'),1]);
%Session = repmat(session,[length(cids'),1]);
unitResponses = table(MouseNum, cids');
unitResponses.Properties.VariableNames = {'MouseNum', 'Cluster'};

% define analysis window
if strcmp(stimuli_parameters.Par.Rec, "SxA")
    PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
    %PostT = 0.25; % captures initial noise period & half of vibrotac (in noise) period
    %PostT = 0.5; % whole vibrotac + dual mode period
    PostT = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
elseif strcmp(stimuli_parameters.Par.Rec, "SOM") && strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'Square')
    PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
    PostT = 0.1; % best way to capture onset stimulus
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

for cond = 1:length(conditions)
    [responsive, hval, pval] = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, conditions(cond));
    
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

    responsive = responsiveCells(stimuli_parameters, baselineRate, stimulusRate, cids, 'OA'); % units responsive

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

%% 3.5 Save if needed
filename = sprintf('M%.2i_S%.2i_%s_ResponseProperties', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
save(fullfile(OutPath, filename), "StimResponseFiring")

%% 4.1 add SOM responses to table
animal = 9; % M11
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
session = [4, 6]; % SxA, SOM M9

OutPath = ['D:\DATA\Processed\M' num2str(animal, '%d') '\ResponseProperties'];

% load SxA data
sessionFile = ['M' num2str(animal, '%.2d') '_S' num2str(session(1), '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
stim_files = dir(fullfile(OutPath, sessionFile));
dataS = load([stim_files.folder '\' stim_files.name]);
unitResponses_SxA = dataS.StimResponseFiring.unitResponses;

% load pressure data
sessionFile = ['M' num2str(animal, '%.2d') '_S' num2str(session(2), '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
stim_files = dir(fullfile(OutPath, sessionFile));
dataS = load([stim_files.folder '\' stim_files.name]);
unitResponses_SOM = dataS.StimResponseFiring.unitResponses;

% combine tables
unitResponses = joindata(unitResponses_SxA, unitResponses_SOM);

%include condition when saving
filename = sprintf('M%.2i_S%.2i&S%.2i_UnitResponses', animal, session(1), session(2));
save(fullfile(OutPath, filename), "unitResponses")

%% 4.2 Bandwith FRA analysis

% parameters
deltaInt = [0,10,20,30,40];
uparamA = unique(stimuli_parameters.Stm.Intensity);
uparamB = unique(stimuli_parameters.Stm.Freq);

% initiate variables
spontRate = StimResponseFiring.firing_mean(1, 1, 1, cluster);
excBand = NaN(length(deltaInt), 2, NClu);
inhBand = NaN(length(deltaInt), 2, NClu);
excBW_oct = NaN(length(deltaInt),NClu);
inhBW_oct = NaN(length(deltaInt),NClu);
excQ = NaN(length(deltaInt),NClu);
inhQ = NaN(length(deltaInt),NClu);
CF = NaN(NClu,1);
BestFreq = NaN(NClu,1);
BestFreqAmp = NaN(NClu,1);

for cluster = 1:NClu

    % find minimum threshold
    [row, col] = find(StimResponseFiring.hvalue(:,:,cluster)==1); % select from sig responses
    if isempty(col)
        minThr = NaN;
    else
        minThr = uparamA(min(col)); % minimal threshold
    end

    %find CF
    CF(cluster) = mean(uparamB(row(col == min(col))));

    % find BF
    [~,maxIdx] = max(StimResponseFiring.firing_mean(:, :, 1, cluster),[],'all','linear'); % find max firing rate change
    [maxRow, maxCol] = ind2sub(size(StimResponseFiring.firing_mean(:, :, 1, cluster)), maxIdx); % translate to position in matrix
    if StimResponseFiring.hvalue(maxRow, maxCol, cluster) % select from sig responses
        BestFreq(cluster) = uparamB(maxRow);
        BestFreqAmp(cluster) = uparamA(maxCol);
    end

    % bandwith
    for dInt = 1:length(deltaInt)
        idx = find(uparamA == minThr + deltaInt(dInt)); % min threshold + delta

        excBandIdx = (StimResponseFiring.firing_mean(idx,:,1,cluster) >= spontRate) & (StimResponseFiring.hvalue(idx,:,cluster)==1);
        inhBandIdx = (StimResponseFiring.firing_mean(idx,:,1,cluster) <= spontRate) & (StimResponseFiring.hvalue(idx,:,cluster)==1);

        if(sum(excBandIdx)>0) % excitatory band
            band_start = find(excBandIdx,1,'first');
            excBand(dInt,1, cluster) = sqrt(uparamB(band_start) * uparamB(band_start - 1)) ;
            band_end = find(excBandIdx,1,'last');
            excBand(dInt,2, cluster) = sqrt(uparamB(band_end) * uparamB(band_end + 1)) ;

            excBW_oct(dInt, cluster) = log2(excBand(dInt,2,cluster) / excBand(dInt,1,cluster));
            excQ(dInt, cluster) = CF(cluster) / (excBand(dInt,2, cluster) - excBand(dInt,1, cluster));
        end

        if(sum(inhBandIdx)>0) % inhibitory band
            band_start = find(inhBandIdx,1,'first');
            inhBand(dInt, 1, cluster) = sqrt(uparamB(band_start) * uparamB(band_start - 1)) ;
            band_end = find(inhBandIdx, 1,'last');
            inhBand(dInt, 2, cluster) = sqrt(uparamB(band_end) * uparamB(band_end + 1)) ;

            inhBW_oct(dInt, cluster) = log2(inhBand(dInt, 2, cluster) / inhBand(dInt, 1, cluster));
            inhQ(dInt, cluster) = CF(cluster) / (inhBand(dInt, 2, cluster) - inhBand(dInt, 1, cluster));
        end
    end
end

%% 4.3 ISI characterization
% for vibrotactile and AM

% inter spike interval analysis
[trials, clusters] = size(aligned_spikes);
aligned_spikes_ISI = cell(trials, clusters);
% set analysis window
if strcmp(stimuli_parameters.Par.Rec, "SxA") && strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
    % vibrotactile stimuli
    stimOn_ISIwindow = str2num(stimuli_parameters.Par.SomAudSOA)/1000;
    stimOff_ISIwindow = stimOn_ISIwindow + str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
elseif  strcmp(stimuli_parameters.Par.Rec, "AMn")
    % AM noise
    stimOn_ISIwindow = ((str2double(stimuli_parameters.Par.AMPostTime)) - (str2num(stimuli_parameters.Par.AMTransTime)))/1000;
    stimOff_ISIwindow = stimOn_ISIwindow + (str2num(stimuli_parameters.Par.AMTransTime))/1000;
end

for cluster = 1:clusters
    for trial = 1:trials
        % select stim on period only from aligned_spikes
        sel = aligned_spikes{trial, cluster} > stimOn_ISIwindow & aligned_spikes{trial, cluster} < stimOff_ISIwindow;
        taligned{trial, cluster} = aligned_spikes{trial,cluster}(sel);

        % log interspike interval
        aligned_spikes_ISI{trial, cluster} = log(diff(taligned{trial, cluster}));
    end
end

clear trial trials cluster clusters sel stimOn_ISIwindow stimOff_ISIwindow taligned

%% 5. combine datasets across mice

% animals to include
% load in stimOrder.csv
stimOrder = readtable('D:\DATA\Processed\stimOrder.csv');

for condition = 2:size(stimOrder, 2)
    
    %initiate variables
    unitResponses_all = table;
    data_all = [];
    mousenum_all =[];
    session_all = [];
    amplitudes_all = [];
    frequencies_all = [];
    PreT_all = [];
    PostT_all = [];

    disp(['Concatinating ' stimOrder.Properties.VariableNames{condition}])

    % get data from session to analyse
    for animal = 1:size(stimOrder,1)

        disp(['Mouse ' num2str(stimOrder{animal,1})])


        if isnan(stimOrder{animal,condition})
            continue
        end

        OutPath = ['D:\DATA\Processed\M' num2str(stimOrder{animal,1}, '%d') '\ICX\ResponseProperties\'];
        sessionFile = ['M' num2str(stimOrder{animal,1}, '%.2d') '_S' num2str(stimOrder{animal,condition}, '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
        stim_files = dir(fullfile(OutPath, sessionFile));
        dataS = load([stim_files.folder '\' stim_files.name]);

        %disp(['post stim time: ' num2str(dataS.StimResponseFiring.PostT) 'sec'])
        %disp(['pre stim time: ' num2str(dataS.StimResponseFiring.PreT) 'sec'])

        % concatinate data from all units
        % select common amp and freq data points
        % if strcmp(dataS.StimResponseFiring.type, 'SxA')
        % 
        %     if animal(animal) == 9 && session(animal) == 4
        %         data = dataS.StimResponseFiring.firing_mean(:, [1,3,5], :, :);
        %     elseif animal(animal) == 8 || animal(animal) == 10 || animal(animal) == 11
        %         data = dataS.StimResponseFiring.firing_mean(1:8, :, :, :);
        %     else
        %         data = dataS.StimResponseFiring.firing_mean; % no changes needed
        %     end
        % else
        %     data = dataS.StimResponseFiring.firing_mean; % no changes needed
        % end

        data = dataS.StimResponseFiring.firing_mean;
        data_all = cat(4, data_all, data);

        % combine unitResponses tables
        % if strcmp(dataS.StimResponseFiring.type, 'SxA')
        %     sessionFile = ['M' num2str(animal(file), '%.2d') '_S' num2str(session(file), '%.2d') '*UnitResponses.mat']; % select based on stimulus type '_*.mat'
        %     stim_files = dir(fullfile(OutPath, sessionFile));
        %     unitResponses = load([stim_files.folder '\' stim_files.name]);
        %     unitResponses_all = vertcat(unitResponses_all, unitResponses.unitResponses);
        % else
        unitResponses_all = vertcat(unitResponses_all, dataS.StimResponseFiring.unitResponses);

        % add details
        mousenum_all = cat(2, mousenum_all, str2double(dataS.StimResponseFiring.MouseNum));
        session_all = cat(2, session_all, str2double(dataS.StimResponseFiring.session));
        amplitudes_all = cat(2, amplitudes_all, dataS.StimResponseFiring.amplitudes);
        frequencies_all = cat(2, frequencies_all, dataS.StimResponseFiring.frequencies);
        PreT_all = cat(2, PreT_all, dataS.StimResponseFiring.PreT);
        PostT_all = cat(2, PostT_all, dataS.StimResponseFiring.PostT);

    end

    StimResponseFiring_all.unitResponses = unitResponses_all;
    StimResponseFiring_all.firing_mean = data_all;
    StimResponseFiring_all.MouseNum = mousenum_all;
    StimResponseFiring_all.session = session_all;
    StimResponseFiring_all.amplitudes = amplitudes_all;
    StimResponseFiring_all.frequencies = frequencies_all;
    StimResponseFiring_all.PreT = PreT_all;
    StimResponseFiring_all.PostT = PostT_all;

    % save
    filename = sprintf('M12-16_%s_UnitResponses',stimOrder.Properties.VariableNames{condition});
    save(fullfile('D:\DATA\Processed\', filename), "StimResponseFiring_all")

    clear data_all

end

%% 5.5 save if needed
filename = sprintf('M09-10-11_SxA_earplug_UnitResponses');
save(fullfile(OutPath, filename), "unitResponses_all", "data_all")

%% ------------------------- Local functions ------------------------- %%
function [cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, onsetDelay, StimResponseFiring] = loadData(OutPath, session, BehaviorPath)
% get data files

% load cluster unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']); %.name;
if ~isempty(cpos_file)
    cpos = load([OutPath '\' cpos_file.name]);
    cids = cpos.clusterinfo.id';
else
    cids = [];
    warning('no cluster details found')
end

% load corresponsing  behaviour files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

% load aligned spike times
aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
if ~isempty(aligned_spikes_files)
    aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
    
    if isfield(aligned_spikes,"Srise")
        disp("Srise/Sfall loaded from data file.")
        Srise = aligned_spikes.Srise;
        Sfall = aligned_spikes.Sfall;
    else
        Srise = [];
        Sfall = [];
        warning("Extract Srise before continuing")
    end

    aligned_spikes = aligned_spikes.SpkT;

else
    Srise = [];
    Sfall = [];
    aligned_spikes = [];
    warning('no aligned spikes file found')
end

% load sessions details
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']); %.name;
if ~isempty(TTLs_file)
    sessions_TTLs = load([OutPath '\' TTLs_file.name]);
    sessions_TTLs = sessions_TTLs.sessions_TTLs;
else
    sessions_TTLs = [];
    warning('no session TTL file found')
end

% load response properties if already exists
if isfolder([OutPath '\ResponseProperties'])
    files = ['\*_S' num2str(session, '%.2d') '_*.mat'];
    resp_files = dir(fullfile([OutPath '\ResponseProperties'], files));
    if ~isempty(resp_files)
        StimResponseFiring = load([resp_files.folder '\' resp_files.name]);
        StimResponseFiring = StimResponseFiring.StimResponseFiring;
        disp('analysed data loaded')
    else
        StimResponseFiring = [];
    end
else
    StimResponseFiring = [];
end

% add delay to "SO"
if strcmp(stimuli_parameters.Par.Rec, 'SxA')
    onsetDelay = stimuli_parameters.Stm.SomAudSOA./1000;
    onsetDelay(isnan(stimuli_parameters.Stm.SomAudSOA)) = 0;
else
    onsetDelay = zeros(size(stimuli_parameters.Stm, 1), 1);
end

end
