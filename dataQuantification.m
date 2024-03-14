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
session = 11;

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

% calculate firing rates pre and post simulus onset

% onset stimuli
PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/4)/1000; % baseline period
PostT = 0.05; % post stim period

% calculate baseline firing rate
[baselineRate, stimulusRate] = baselinerate(aligned_spikes.SpkT, PreT, PostT); % get firing rates (Hz) in time window

%% plot summary fig
% boxplot, per unit, vibrotac sessions

close all
amp = 0.3; 

uFreq = unique(stimuli_parameters.Stm.SomFreq);
nFreq = length(uFreq);
nClusters = length(cids);
condition = nan(nClusters, nFreq);
figure

for cluster = 1:nClusters
    for freq = 1:nFreq
        condition(cluster, freq) = median(stimulusRate([stimuli_parameters.Stm.SomFreq] == uFreq(freq) & [stimuli_parameters.Stm.Amplitude] == amp, cluster))';
    end

    plot(uFreq, mean(condition(cluster, :) - condition(cluster, 1), 1), ':o')
    hold on

end

xticks(uFreq)
xticklabels(uFreq(1:nFreq));
xlabel(['Vibrotactile stimulus (Hz), ' num2str(amp) ' (V)'])
ylabel('Normalized firing rate (Hz)')

title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

save figure
figname = sprintf('M%.2i_S%.2i_%s_mean', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

% plot summary fig: boxplot, all units, vibrotac sessions
close all
condition = nan(nClusters, nFreq);
figure

% select data
for freq = 1:nFreq
    condition(:, freq) = median(stimulusRate([stimuli_parameters.Stm.SomFreq] == uFreq(freq) & [stimuli_parameters.Stm.Amplitude] == amp, :))';
end

% plot
figure
boxplot(condition)

xticklabels(uFreq(1:nFreq));
xlabel(['Vibrotactile stimulus (Hz), ' num2str(amp) ' (V)'])
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

% save figure
figname = sprintf('M%.2i_S%.2i_%s_boxplot_all units', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));


%% plot summary fig
% boxplot, per  units, pressure sessions
close all

% parameters
uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
nClusters = length(cids);
condition = nan(nClusters, nAmp);

figure
for cluster = 1:nClusters

    % select data
    for amplitude = 1:nAmp
        condition(cluster, amplitude) = median(stimulusRate([stimuli_parameters.Stm.Amplitude] == uAmp(amplitude), cluster))';
    end

    % plot
    plot(uAmp, mean(condition(cluster, :) - condition(cluster, 1), 1), ':o');
    hold on

end

xticks(uAmp)
xticklabels(uAmp(1:nAmp));
xlabel('Pressure (V)')
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

%save figure
figname = sprintf('M%.2i_S%.2i_%s_mean', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

% plot summary fig: boxplot, all units, pressure sessions
close all
condition = nan(nClusters, nAmp);

% select data
for amplitude = 1:nAmp
    condition(:, amplitude) = median(stimulusRate([stimuli_parameters.Stm.Amplitude] == uAmp(amplitude), :))';
end

% plot
figure
boxplot(condition)

xticklabels(uAmp(1:nAmp));
xlabel('Pressure (V)')
ylabel('Normalized firing rate (Hz)')
title(['Stimulus induced responses - session ' stimuli_parameters.Par.Set ': ' stimuli_parameters.Par.SomatosensoryLocation])

%save figure
figname = sprintf('M%.2i_S%.2i_%s_boxplot_all units', str2double(stimuli_parameters.Par.MouseNum), str2double(stimuli_parameters.Par.Set), stimuli_parameters.Par.Rec);
saveas(gcf, fullfile(OutPath, [figname '.jpg']));
saveas(gcf, fullfile(OutPath, figname));

%% plot firing rate differences

% select correct part of matrix
% off all units: median firing rate after stim
% iqr

   
%% quantify reactive units
% to do
% 1. quantify responsive units: sig diff firing rate between stim period vs baseline

% onset stimuli
PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/5)/1000; % baseline period
PostT = 0.05; % post stim period

uAmp = unique(stimuli_parameters.Stm.Amplitude);
nAmp = length(uAmp);
uFreq = unique(stimuli_parameters.Stm.SomFreq);
nFreq = length(uFreq);
nClusters = length(cids);

results.cids = cids;
results.conditions = uAmp;
results.pvalue = nan(nAmp, nClusters);
results.signrank = nan(nAmp, nClusters);
results.hvalue = nan(nAmp, nClusters);
results.responsive = [];

% stats test difference baseline vs stimulus period
for cluster = 1:nClusters

    % calculate baseline firing rate
    [baselineRate, stimulusRate] = baselinerate(aligned_spikes.SpkT, PreT, PostT); % get firing rates (Hz) in time window
   
    %signrank(baselineRate(:,cluster), stimulusRate(:,cluster), 'alpha', 0.01); % overall different from baseline?
    control = stimuli_parameters.Stm.Amplitude == 0;

    for condition = 1:nAmp
        index = stimuli_parameters.Stm.Amplitude == uAmp(condition);
        [p,h,stats] = signrank(stimulusRate(control,cluster), stimulusRate(index,cluster), 'alpha', 0.01); % different from control trials?
        results.pvalue(condition, cluster) = p;
        results.signrank(condition, cluster) = stats.signedrank;
        results.hvalue(condition, cluster) = h;
    end

    % unit responsive if sig diff for at least one condition
    if max(results.hvalue(:, cluster)) == 1
        results.responsive = [results.responsive, results.cids(cluster)];
    end

end


%% stim session
% 2. cross correlating single trials (KDF as in previous script), corr
% coeficient over control value to be sig
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
        
        % get spike times per trial
        S = SpkT{trial, cluster};
        
        % count spikes in time windows
        baseSpikes = sum(S <= 0);
        expSpikes = sum(S > 0 & S <= PostT);

        % convert to firing rate (Hz)
        baseFR(trial, cluster) = baseSpikes/abs(PreT);
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
