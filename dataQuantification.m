% To be done after PostCuration_BK.m
% General data input: single units (array) stimuli order (table), aligned spikes (cell array),
% Post processing data steps:
%   1. get firing rate during set time frame for each trial of each unit
%   2. quatify which cluster responds to which event
%       responsive unit = sig different firing after stim vs baseline
%   3. quantify response strength
%   4. combine data from multiple animals

clearvars

%SxA, SOM session matching 
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
%session = [4, 6]; % M9

%set paths
mousenr = 19;
OutPath = 'D:\DATA\Processed\M19\ICX';
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M19';

%% quantitative analysis

sessions = [1 2 6 7 8 9 10 11 12 13 14 15 16]; %all sessions to be analyzed 
for session = 1:length(sessions)
    if isempty(dir([OutPath '\*_S' num2str(sessions(session),'%02i') '*_ResponseProperties.mat']))
        dataQuantification_analysis(mousenr, OutPath, BehaviorPath, sessions(session))
    else
        continue
    end
end
%% Pool datasets across mice
stimOrder = readtable('D:\DATA\Processed\M10-11-19-20_stimOrder.csv');
fn = 'M10-11-19-20';
%%
dataQuantification_poolDatasets(stimOrder, fn)

%% add mousenr to unit nr

Path = 'D:\DATA\Processed\M10-11-19-20\';
File = '*_UnitResponses.mat';
files = dir(fullfile(Path, File));

for i = 3:size(files, 1)
    data = load([files(i).folder '\' files(i).name]);

    cids = data.StimResponseFiring_all.unitResponses.Cluster';
    mice = data.StimResponseFiring_all.unitResponses.MouseNum';
    s = sprintf( '%1d%1d ', [mice;cids]);
    data.StimResponseFiring_all.unitResponses.Cluster = str2num(s)';

    StimResponseFiring_all = data.StimResponseFiring_all;

    %filename = sprintf('M10-11-19-20_%s_UnitResponses', stimOrder.Properties.VariableNames{condition});
    save(fullfile(files(i).folder, files(i).name), "StimResponseFiring_all")
end

%% STATS - to do

%% OPTIONAL: add SOM responses to SxA table, if needed
mousenr = 9; % M11
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
session = [4, 6]; % SxA, SOM M9

OutPath = ['D:\DATA\Processed\M' num2str(mousenr, '%d') '\ResponseProperties'];

% load SxA data
sessionFile = ['M' num2str(mousenr, '%.2d') '_S' num2str(session(1), '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
stim_files = dir(fullfile(OutPath, sessionFile));
dataS = load([stim_files.folder '\' stim_files.name]);
unitResponses_SxA = dataS.StimResponseFiring.unitResponses;

% load pressure data
sessionFile = ['M' num2str(mousenr, '%.2d') '_S' num2str(session(2), '%.2d') '_*_ResponseProperties.mat']; % select based on stimulus type '_*.mat'
stim_files = dir(fullfile(OutPath, sessionFile));
dataS = load([stim_files.folder '\' stim_files.name]);
unitResponses_SOM = dataS.StimResponseFiring.unitResponses;

% combine tables
unitResponses = joindata(unitResponses_SxA, unitResponses_SOM);

%include condition when saving
filename = sprintf('M%.2i_S%.2i&S%.2i_UnitResponses', mousenr, session(1), session(2));
save(fullfile(OutPath, filename), "unitResponses")

%% OPTIONAL: Bandwith FRA analysis

session = 4;
[cids, stimuli_parameters, ~, ~, ~, ~, ~, StimResponseFiring] = loadData(OutPath, session, BehaviorPath);

% parameters
deltaInt = [0,10,20,30,40];
uparamA = unique(stimuli_parameters.Stm.Intensity);
uparamB = unique(stimuli_parameters.Stm.Freq);
NClu = length(cids);

% initiate variables
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
    spontRate = StimResponseFiring.firing_mean(1, 1, 1, cluster);

    % find minimum threshold
    [row, col] = find(StimResponseFiring.hvalue(:,:,cluster)==1); % select from sig responses
    if isempty(col)
        minThr = NaN;
    else
        minThr = uparamA(min(col)); % minimal threshold
    end

    %find charactristic frequency (CF)
    CF(cluster) = mean(uparamB(row(col == min(col))));

    % find best frequency (BF)
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
[cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, StimResponseFiring] = loadData(OutPath, 3, BehaviorPath);

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

conditions = StimResponseFiring.conditions; % OO OA SO SA
uparamA = StimResponseFiring.amplitudes;
uparamB = StimResponseFiring.frequencies;
nparamA = length(uparamA);
nparamB = length(uparamB);
ISICVs = nan(nparamA, nparamB, length(conditions), length(cids)); % most freq ISI %amp,freq,con,cluster

spikeISI = cell(nparamA, nparamB, size(conditions, 2), length(cids));
for cluster = 1:length(StimResponseFiring.cids)
    for con = 1:length(conditions)
        for amp = 1:nparamA
            for freq = 1:nparamB

                % select unique stimulus combination
                if strcmp(stimuli_parameters.Par.Rec, "SxA")
                    index = stimuli_parameters.Stm.Amplitude == uparamA(amp) & stimuli_parameters.Stm.SomFreq == uparamB(freq) & strcmp(stimuli_parameters.Stm.MMType, conditions{con});
                    figTitle = [num2str(uparamB(freq)) 'Hz ' num2str(uparamA(amp)) 'V'];
                    figSGTitle = ['Unit: ' num2str(cids(cluster)) ', condition: ' conditions{con} '(' num2str(uparamA(amp)) 'V)'];
                elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
                    index = stimuli_parameters.Stm.Intensity == uparamA(amp) & stimuli_parameters.Stm.Mf == uparamB(freq);
                    figTitle = [num2str(uparamB(freq)) 'Hz '];
                    figSGTitle = ['Unit: ' num2str(cids(cluster)) ', ' num2str(uparamA(amp)) 'dbSPL'];
                end

                if sum(index) == 0
                    continue
                end

                idx_rows = find(index);
                tspikeISI = [];

                % concatinate all ISIs
                for i = 1:sum(index)
                    ttISIspikes = aligned_spikes_ISI{idx_rows(i), cluster};
                    tspikeISI = [tspikeISI; ttISIspikes];
                end

                % to do: save SpikeISI in cell array to disentangle analysis from plotting
                spikeISI{amp, freq, con, cluster} = tspikeISI; %amp, freq, con, cluster

                % calculate coefficient of variation (CV = standard deviation / mean)
                ISICV = mean(exp(spikeISI{amp, freq, con, cluster}))/ std(exp(spikeISI{amp, freq, con, cluster}));
                ISICVs(amp,freq,con,cluster) = ISICV;

                % find most occuring ISI
                ISImax = max(N);
                ISImaxpos = find(N == ISImax);

                if length(ISImaxpos) == 1
                    ISIpeak = edges(ISImaxpos); % match with edges to get ISI
                else
                    %disp('two peaks found, check')
                    ISIpeak = edges(ISImaxpos(1)); % match with edges to get ISI
                end

                % save in matrix amp x freq x condition x cluster
                ISIpeaks(amp,freq,con,cluster) = ISIpeak;

            end
        end
    end
end





clear trial trials cluster clusters sel stimOn_ISIwindow stimOff_ISIwindow taligned tspikeISI ttspikeISI
