%% Quantify results
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   quatify which cluster responds to which event
%       response = fire (spiketime) in event window
%   quantify firing rates? changes expected after events vs spontaneous

clearvars

% set directories
%recordingFolder = 'D:\DATA\EphysRecordings\M8\M08_2024-02-27_12-29-52\Record Node 103\experiment1\recording1\';
%recPath = [recordingFolder 'continuous\Intan-100.Rhythm Data-A\'];
%TTLPath = [recordingFolder 'events\Intan-100.Rhythm Data-A\TTL\'];
%messagesPath = [recordingFolder 'events\MessageCenter\'];
%rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording
%KSPath = 'D:\DATA\EphysRecordingsSorted\M08\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M8\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M8'; % output directory

% select which session to analyse
session = 2;

% load unit info
cpos_file = dir([OutPath '\*_InfoGoodUnits.mat']).name;
cpos = load([OutPath '\' cpos_file]);
cids = cpos.cpos.id';

% load sessions details
TTLs_file = dir([OutPath '\*_OE_TTLs.mat']).name;
sessions_TTLs = load([OutPath '\' TTLs_file]);

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile(OutPath, sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);

%% quantify reactive units
% to do
% 1. quantify responsive units: sig diff firing rate between stim period vs baseline

% calculate baseline firing rate
PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/4)/1000; % baseline period
PostT = 0.5; % post stim period

[baselineRate, stimulusRate] = baselinerate(aligned_spikes.SpkT, PreT, PostT); % get firing rates (Hz) in time window

uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
nFreq = unique(stimuli_parameters.Stm.SomFreq);

p = nan(condition, cluster);
h = nan(condition, cluster);
stats = nan(condition, cluster);

for cluster = 1:length(cids)
    signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    control = stimuli_parameters.Stm.Amplitude == 0;

    for condition = 1:nAmp
        index = stimuli_parameters.Stm.Amplitude == uAmp(condition);
        [p,h,stats] = signrank(stimulusRate(control,cluster), stimulusRate(index,cluster), 'alpha', 0.01); % different from control trials?
        p(condition, cluster) = p;
        h(condition, cluster) = h;
        stats(condition, cluster) = stats;
    end
end

%% stim session
% 2. cross correlating single trials (KDF as in previous script), corr
% coeficient over control value to be sig
% onset/offset responses quantify with window (0.1sec?)
% 3. firing in relation to phase stimulus (try cycle histogram)

%%
SponR = spontaneousrate(SpkS, S_sets);


%% calculate baseline firing rate
function [baseFR, expFR] = baselinerate(SpkT, PreT, PostT)

nClust = size(SpkT,2);
nTrials = size(SpkT,1);
baseFR = nan(nTrials, nClust);
expFR = nan(nTrials, nClust);

for cluster = 1:size(SpkT,2)
    for trial=1:size(SpkT,1)
        S = SpkT{trial, cluster};
        baseSpikes = sum(S <= 0);
        expSpikes = sum(S > 0 & S <= PostT);
        baseFR(trial, cluster) = baseSpikes/abs(PreT); % count spikes before stim on
        expFR(trial, cluster) = expSpikes/abs(PostT);
    end
end
end

%%

function [sponR] = spontaneousrate(SpkS,S_sets)

USets = S_sets.setNum;
nSets = S_sets.NSets;
paths = S_sets.path;

nClust = size(SpkS,2);

sponR = nan(nSets,nClust);

% Srise = Stm.Srise;

for k=1:nSets
    File = paths{1,k};
    Rec = string(File(1,end-6:end-4));
    if Rec == "FRA" || Rec == "AMn"
        load(paths{1,k},'Stm');
        Srise = Stm.Srise;
        Srise = Srise(1);
    elseif Rec == "DRC" || Rec == "OMI"
        load(paths{1,k},'DigSamp');
        Srise = DigSamp{1,2}; Srise=Srise(1,1);
    end
    
  
    for s=1:nClust
        S = SpkS{USets(k),s};
        sponSpkS = S(S<=Srise); %extract all spikes prior to the start of the first trial of the set
        sponR(k,s) = (length(sponSpkS)/(Srise/30000)); %divide the spike count from this period by measurement period and store in sepSFR
        
    end
    
end

end
