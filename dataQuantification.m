% To be done after PostCuration_BK.m
% General data input: single units (array) stimuli order (table), aligned spikes (cell array),
% Post processing data steps:
%   1. get firing rate during set time frame for each trial of each unit
%   2. quatify which cluster responds to which event
%       responsive unit = sig different firing after stim vs baseline
%   3. quantify response strength
%   4. combine data from multiple animals

%SxA, SOM session matching 
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
%session = [4, 6]; % M9

% set paths
%OutPath
%session
%BehaviorPath


%% quantitative analysis
dataQuantification_analysis(animal, OutPath, BehaviorPath, session)
%filename = sprintf('M%.2i_S%.2i_%s_ResponseProperties', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
%save(fullfile(OutPath, filename), "StimResponseFiring")

%% Pool datasets across mice
stimOrder = readtable('D:\DATA\Processed\stimOrder.csv');
dataQuantification_poolDatasets(stimOrder, fn)

%% STATS - to do

%% OPTIONAL: add SOM responses to SxA table, if needed
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

%% OPTIONAL: Bandwith FRA analysis

[~, stimuli_parameters, ~, ~, ~, ~, ~, StimResponseFiring] = loadData(OutPath, session, BehaviorPath);

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

%% OPTIONAL: ISI characterization
% for vibrotactile and AM

% load relevant data
[~, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

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


