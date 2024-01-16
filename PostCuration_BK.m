%% Post-processing ephys data
% data paths to sorted data, behavioral data file, output folder
% post processing data steps:
%   get spike times from all good clusters
%   get event timings (from ephys events)
%   quatify which cluster responds to which event
%       response = fire (spiketime) in event window
%   quantify firing rates? changes expected after events vs spontaneous

clearvars

% set directories
recPath = 'D:\DATA\EphysRecordings\M4\M04_2023-12-14_12-47-41\Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data\';
TTLPath = 'D:\DATA\EphysRecordings\M4\M04_2023-12-14_12-47-41\Record Node 103\experiment1\recording1\events\Intan-100.Rhythm Data\TTL\';
messagesPath = 'D:\DATA\EphysRecordings\M4\M04_2023-12-14_12-47-41\Record Node 103\experiment1\recording1\events\MessageCenter\'; % session TTLs
KSPath = 'D:\DATA\EphysRecordingsSorted\M04\trimmed\'; % kilosort ephys data
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M4\'; % stimuli parameters
OutPath = 'D:\DATA\Processed\M4'; % output directory

rec_samples = readNPY([recPath 'sample_numbers.npy']); % sample nr whole recording

relevant_sessions = [1 23]; % behaviour files (if only 1 behavior file in rec: [1 1])

%% IronClust: post-curation unit extraction

%[spiketimes, cids,cpos] = ircGoodClusters(spiketimecsv,clusterqualitycsv);

%% Kilosort: post-curation unit extraction
% spike extraction from curated units

[spiketimes, cids, cpos, Srise, Sfall] = extractspikes(BehaviorPath, KSPath, recPath, TTLPath, messagesPath, relevant_sessions, rec_samples);

%% save details good units
set = sprintf('%02d-%02d', relevant_sessions(1), relevant_sessions(2));
filename = ['M04_S' set '_InfoGoodUnits'];
save(fullfile('D:\DATA\Processed', filename), "cpos") %cpos variables: unit id, channel, depth, avg firing rate, nr spikes;

%% align spikes
% alignment of extracted spikes to stimulus on/off-set
% spike times in sec
% good to have: TTL signaling start & end of stim presentation set

% function saves aligned spikes cell array, change mouse name
alignspikes(BehaviorPath, spiketimes, relevant_sessions, Srise, Sfall, cids);

%% FRA analysis
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster

% select correct input files
aligned_spikes = load('D:\DATA\Processed\M04_S01_FRA_AlignedSpikes');
stimuli_parameters = load([BehaviorPath 'M4_S01_FRA.mat']);

% function saves figures, change mouse name
FRAanalysis(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath, 1);

%% plotting single sessions
% FRA: raster
% AMn: raster + PSTH
% SOM: raster + PSTH

%cids = load('D:\DATA\Processed\M04_S01-23_InfoGoodUnits.mat'); 

% select which session to plot
session = 22;

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile('D:\DATA\Processed', sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);

plotResponses(stimuli_parameters, aligned_spikes.SpkT, cids, OutPath);

%% plotting SOM sessions with their controls
%saveplots = 0; %0 don't save, 1 save plots in OutPath
sessions = [23 21]; % [exp ctrl]
%cids = load('D:\DATA\Processed\M04_S01-23_InfoGoodUnits.mat'); 
cids = [    90   111   124   126   151   159   169   171   182   188   218   232   256   261   264   265   268   278];
OutPath = 'D:\DATA\Processed\M4';

SOMplotting(sessions, cids, OutPath, BehaviorPath, 0);
%% FSL SOM analysis
% function SOM = SOManalysis(stimuli_parameters, aligned_spikes, cids)
% input: stimuli_parameters.Par, stimuli_parameters.Stm, aligned_spikes
% output: first spike latency SOM/AM trials

% select which session to plot
session = 23;

% load corresponsing files
sessionFile = ['\*_S' num2str(session, '%.2d') '_*.mat'];
stim_files = dir(fullfile(BehaviorPath, sessionFile));
stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

aligned_spikes_files = dir(fullfile('D:\DATA\Processed', sessionFile));
aligned_spikes = load([aligned_spikes_files.folder '\' aligned_spikes_files.name]);
aligned_spikes = aligned_spikes.SpkT;

cids = [90 111 124 126 151 159 169 171 182 188 218 232 256 261 264 265 268 278];

SOM = FSL_SOM_AMn(stimuli_parameters, aligned_spikes, cids);

%% plotting channel map
% match unit position to channel map

channel_map = readNPY([KSPath 'channel_map.npy']);
channel_positions = readNPY([KSPath 'channel_positions.npy']);
cids = load('D:\DATA\Processed\M04_S01-23_InfoGoodUnits.mat'); 
cpos = cids.cpos;

xcoords = channel_positions(:,1);
ycoords = channel_positions(:,2);
margin = 50;
xMin = min(xcoords); xMax = max(xcoords);
xSpan = xMax-xMin; xMid = 0.5*(xMax+xMin);
yMin = min(ycoords); yMax = max(ycoords);
ySpan = yMax-yMin; yMid = 0.5*(yMax+yMin);
maxSpan = max(xSpan,ySpan);

fig = scatter(xcoords, ycoords, ".", 'k');
hold on;

for i = 1:length(channel_map)
    text((channel_positions(i,1)+1), (channel_positions(i,2)+1),num2str(channel_map(i)))
end

% add position of units included in analysis
hold on
units = table2array(cpos(:,2));
idx = ismember(channel_map, units);
fig = scatter(channel_positions(idx,1), channel_positions(idx,2), 'o'); % wrongly indexed

% to do: reflect number of clusters on 1 channel. apply jitter in circle
% marking to achieve this
% to do: add unit number with corresponding channel
%text((channel_positions(unique(idx),1)+1), (channel_positions(unique(idx),1)+1),num2str(channel_map()))

% format figure
title('Channel map')
ylabel('Relative depth')
fig.SizeData = 100;
axis square
xlim([xMid - 0.5*maxSpan - margin, xMid + 0.5*maxSpan + margin]);
ylim([yMid - 0.5*maxSpan - margin, yMid + 0.5*maxSpan + margin]);
xticks(unique(xcoords));
yticks(unique(ycoords));


%% channel waveforms
% extract and plot waveform traces
addpath('C:\Users\TDT\Documents\GitHub\ephys_analysis');
for k=1:length(cids)
    close all
    %F = myfig(1,'fig');
    F = figure;
    F.Position = [2405,214,838,890];
    F = ChannelWaveforms(F,KSPath,cids(k),cpos(k,3),200);
    F.Renderer = 'painter';
    F.PaperUnits = 'inches';
    F.PaperSize = F.Position(3:4)./96;
    %saveas(F,[UserPath 'Processed\' DirName '\Units\' DirName '-U' num2str(cids(k)) '.pdf']);
    %     saveas(F,[UserPath 'Processed\' Mouse '\Units\' Mouse '-U' num2str(cids(k)) '.png']);
end

%% quantify reactive units
% to do

%% FRA analysis old
% works correctly if whole recording corresponds with 1 behaviour file

% quantify firing rate per unit in response to every sound intensity x frequency combination
% needed for this:
% aligned_spikes (trial x unit, # spike times per combo) amount of spikes in certain time window
% freq = nrspikes / 290 (ms, timewindow alignment)

% TO DO: select correct part of aligned_spikes to plot
% TO DO: make into function

NClu = length(cids); % cluster info

NSets = 1;
setNum = 1; % set info
%Stm = S_sets.Stm(StmIdx, :); % what does this do?

% stimulus parameters - frecuency
UFreq = unique([stimuli_parameters.Stm.Freq]);
NFreq = length(UFreq);

% stimulus parameters - intensity
UInt = unique([stimuli_parameters.Stm.Intensity]);
UInt = UInt(~isnan(UInt)); % get no NaNs
NInt = length(UInt);

% analysis period
winStart = 0;
winEnd   = 0.115;

% % initiation variables
FRAScntSD = nan(NInt, NFreq, NSets, NClu); % spike count
FRAScnt = nan(NInt, NFreq, NSets, NClu); % spike count
FRASR = nan(NInt, NFreq, NSets, NClu); % spike count
NTrials = nan(NInt, NFreq, NSets, NClu); % number of trials
MedFSL = nan(NInt, NFreq, NSets, NClu); % median first spike latency
% % FACA = nan(NInt, NFreq, NSets, NClu); % correlation
% % pval = nan(NInt, NFreq, NSets, NClu); % correlation
% % AvgNeighCorr = nan(NInt, NFreq, NSets, NClu); % correlation
% % MaxNeighCorr = nan(NInt, NFreq, NSets, NClu); % correlation
% % avgKSD = cell(NInt, NFreq, NSets, NClu); % spike density
% % KSDtt = cell(NSets, 1); % spike density


for cluster = 1:NClu
    disp(['Analysing cluster ' num2str(cluster) ' of ' num2str(NClu)]);

    % Spike count analysis
    for freq = 1:NFreq
        for intensity = 1:NInt
            sel = [stimuli_parameters.Stm.Freq] == UFreq(freq) & [stimuli_parameters.Stm.Intensity] == UInt(intensity);
            NTrials(intensity,freq,setNum,cluster) = sum(sel);

            if (NTrials(intensity,freq,setNum,cluster) == 0); continue; end

            %count spikes
            SCnt = nan(NTrials(intensity,freq,setNum,cluster), 1);

            tempSpiketimes = aligned_spikes(sel,cluster);
            for t = 1:length(tempSpiketimes)
                if (isnan(tempSpiketimes{t})); continue; end
                if(isempty(tempSpiketimes{t}))
                    SCnt(t) = 0;
                else
                    SCnt(t) = sum(tempSpiketimes{t} > winStart & tempSpiketimes{t} < winEnd);
                end
            end

            %SCnts = countSpks(spiketimes(sel), winStart, winEnd);
            %SCnts = countSpks(tempSpkT(sel),tempWinStart,tempWinEnd);
            %SpkT, what does this look like?

            FRAScntSD(intensity,freq,setNum,cluster)  = std(SCnt,'omitnan');
            FRAScnt(intensity,freq,setNum,cluster)    = mean(SCnt,'omitnan');
            FRASR(intensity,freq,setNum,cluster)      = mean(SCnt,'omitnan') / (winEnd - winStart);

            % to do: adapt for SOM
            %FSL function
            fsl = inf(length(tempSpiketimes), 1);
            for t = 1:length(tempSpiketimes)
                if (isnan(tempSpiketimes{t}))
                    continue
                end

                spks = tempSpiketimes{t};
                spks = spks (spks > 0);

                if (~isempty(spks))
                    fsl(t) = min(spks);
                end
            end

            MedFSL(intensity, freq, setNum, cluster)     = median(fsl);

        end % intensity loop
    end % frequency loop
    clearvars('tempSpiketimes');

    % %% Correlation based analysis
    %
    % % Calculating Spike Density
    % NStim = length(tempSpk);
    % tStep = 0.5e-3; %s
    % minT = -min([tempStm.PreT]) * 1e-3; % ms -> s
    % maxT = min([tempStm.StimT]+[tempStm.PostT]) * 1e-3; % ms -> s
    % xx = minT+bw*2:tStep:maxT-bw*2;
    % KSD = nan(length(xx),NStim);
    % for ss = 1:NStim
    %     if isempty(tempSpk{ss})
    %         KSD(:,ss) = 0; continue
    %     end
    %     [KSD(:,ss),~] = ksdensity(tempSpk{ss},xx,'Bandwidth',bw);
    %     KSD(:,ss) = numel(tempSpk{ss}) * KSD(:,ss); % normalize to spike density, i.e. cumsum * dT = Nspikes
    % end
    %
    %
    %
    %
    % % Averaging Spike Density
    % for i = 1:NInt
    %     Int = UInt(i);
    %     for f = 1:NFreq
    %         Freq = UFreq(f);
    %         idx = find([tempStm.Freq] == Freq &  [tempStm.Intensity] == Int);
    %         if isempty(idx); continue; end
    %         avgKSD{i,f,s,c} =  mean(KSD(:,idx),2);
    %     end
    % end
    %
    %
    %
    % % Calculating Correlation Coefficients
    % NStim           = size(KSD,2);
    % rho             = corr(KSD);
    % rho(isnan(rho)) = 0;
    % rho             = rho+diag(nan(NStim,1));
    %
    %
    % % Sampling FACA
    % for i = 1:NInt
    %     Int = UInt(i);
    %     for f = 1:NFreq
    %         Freq = UFreq(f);
    %         idx = find([tempStm.Freq] == Freq &  [tempStm.Intensity] == Int);
    %         if isempty(idx); continue; end
    %         FACA(i,f,s,c) =  mean(rho(idx,idx),'all','omitnan');
    %     end
    % end
    %
    % setUFreq    = unique([tempStm.Freq]);
    % fIdx        = ismember(UFreq,setUFreq);
    % setUInt     = unique([tempStm.Intensity]);
    % iIdx        = ismember(UInt,setUInt);
    %
    % % Sampling NeighbourCorrelation
    % diagWeight = 0.5;
    %
    % UInt = unique([tempStm.Intensity]);
    % NInt = length(UInt);
    % UFreq = unique([tempStm.Freq]);
    % NFreq = length(UFreq);
    %
    % vCorr   = nan(NInt,NFreq);
    % hCorr   = nan(NInt,NFreq);
    % d1Corr  = nan(NInt,NFreq);
    % d2Corr = nan(NInt,NFreq);
    %
    % for ii = 1:NInt
    %     Int = UInt(ii);
    %     for ff = 1:NFreq
    %         Freq = UFreq(ff);
    %         if Freq == 0 || isinf(Int); continue; end % silent trials
    %         sel = [tempStm.Freq] == Freq &  [tempStm.Intensity] == Int;
    %         idx = find(sel);
    %         if isempty(idx); continue; end % empty
    %         % vertical
    %         if ii < NInt
    %             Int2 = UInt(ii+1);
    %             Freq2   = Freq;
    %             sel2 = [tempStm.Freq] == Freq2 &  [tempStm.Intensity] == Int2;
    %             idx2 = find(sel2);
    %             if ~isempty(idx2)
    %                 vCorr(ii,ff) =  mean(rho(idx,idx2),'all','omitnan');
    %             end
    %         end
    %         % horizontal
    %         if ff < NFreq
    %             Int2    = Int;
    %             Freq2   = UFreq(ff+1);
    %             sel2    = [tempStm.Freq] == Freq2 &  [tempStm.Intensity] == Int2;
    %             idx2    = find(sel2);
    %             if ~isempty(idx2)
    %                 hCorr(ii,ff) =  mean(rho(idx,idx2),'all','omitnan');
    %             end
    %         end
    %         % diag1
    %         if ff < NFreq && ii < NInt
    %             Int2    = UInt(ii+1);
    %             Freq2   = UFreq(ff+1);
    %             sel2    = [tempStm.Freq] == Freq2 &  [tempStm.Intensity] == Int2;
    %             idx2    = find(sel2);
    %             if ~isempty(idx2)
    %                 d1Corr(ii,ff) =  mean(rho(idx,idx2),'all','omitnan');
    %             end
    %         end
    %         % diag2
    %         if ff > 1 && ii < NInt
    %             Int2    = UInt(ii+1);
    %             Freq2   = UFreq(ff-1);
    %             sel2    = [tempStm.Freq] == Freq2 &  [tempStm.Intensity] == Int2;
    %             idx2    = find(sel2);
    %             if ~isempty(idx2)
    %                 d2Corr(ii,ff) =  mean(rho(idx,idx2),'all','omitnan');
    %             end
    %         end
    %
    %     end
    % end
    %
    % dU  = vCorr; dD = circshift(vCorr,1,1);
    % dR  = hCorr; dL = circshift(hCorr,1,2);
    % dTL = d1Corr; dBR = circshift(circshift(d1Corr,1,2),1,1);
    % dTR = d2Corr; dBL = circshift(d2Corr,1,1);
    %
    % NeighCorr = cat(3,dU,dD,dR,dL,diagWeight.*dTL,diagWeight.*dBR,diagWeight.*dTR,diagWeight*dBL);
    % AvgNeighCorr = mean(NeighCorr,3,'omitnan');
    % MaxNeighCorr = max(NeighCorr,[],3,'omitnan');
    %
    % % Sampling Average Correlation Coefficient
    %
    % NStim = size(rho,1);
    % NSamp = 50000;
    % if (max(rho,[],'all') == 0); M_Bootstrap = zeros(NSamp,maxRep); return;end
    % M_Bootstrap = nan(NSamp,maxRep);
    % % look up
    % for samp = 1:NSamp
    %     permList = randperm(NStim,maxRep);
    %     for sampSize = 1:maxRep
    %         M_Bootstrap(samp,sampSize) = mean(rho(permList(1:sampSize),permList(1:sampSize)),'all','omitnan');
    %     end
    % end
    % M_Bootstrap = sort(M_Bootstrap);
    %
    %
    % % Looking up p-values
    %
    %
    % alpha = 0.01;
    %
    % pval = nan(size(FACA));
    % UNTrials = unique(NTrials);
    % UNTrials = UNTrials(UNTrials>1);
    % for nn = 1:length(UNTrials)
    %     sel =  NTrials(:) == UNTrials(nn) ;
    %
    %     NSamp = size(M_Bootstrap,1);
    %     [CDF_x,CDF_y,~] = unique(M_Bootstrap(:,SampSize),'last');
    %     CDF_y = (CDF_y - 1) ./ NSamp;
    %     if (length(CDF_x) < 2); pval = nan(size(FACA)); return;end
    %     pval = 1-interp1(CDF_x,CDF_y,FACA);
    %     pval(FACA > max(CDF_x)) = 1/NSamp/2;
    %     pval(FACA < min(CDF_x)) = 1-1/NSamp/2;
    %
    %
    % end
    %
    % sig = double(pval < alpha/2 | pval > 1-alpha/2); % two tail
    % sig(isnan(FACA)) = nan;



end % cluster loop

% Output data
% Experiment meta data
% FRA.FRASets      = FRASets;
% FRA.FRAIdx       = FRAIdx;
% FRA.FRASetNum    = FRASetNum;
FRA.NSets        = NSets;
% FRA.startTime    = S_sets.startTime(FRAIdx);
% FRA.meanTime     = S_sets.meanTime(FRAIdx);

% cluster
FRA.NClu = NClu;

% stimulus parameters
FRA.UFreq = UFreq;
FRA.NFreq = NFreq;
FRA.UInt = UInt;
FRA.NInt = NInt;

% results - number of trials
FRA.NTrials = NTrials;

% results - spike latency
FRA.MedFSL = MedFSL;

% results - spike count
FRA.FRASR = FRASR;
FRA.FRAScnt = FRAScnt;
FRA.FRAScntSD = FRAScntSD;

% results - spike density
% FRA.avgKSD  = avgKSD;
% FRA.KSDtt   = KSDtt;

% results - correlation
% FRA.FACA        = FACA;
% FRA.FACApval    = pval;
% % FRA.FACAsig     = sig;
% FRA.AvgNeighCorr = AvgNeighCorr;
% FRA.MaxNeighCorr = MaxNeighCorr;


%% FRA heatmap
for clustNum = 1:NClu
    NSets       =	1;
    NClu       =	FRA.NClu;
    UFreq       =	FRA.UFreq;
    NFreq       =	FRA.NFreq;
    UInt       =	FRA.UInt;
    NInt       =	FRA.NInt;

    startTime = 0;%FRA.startTime;
    meanTime = 0;%FRA.meanTime;

    M   =   FRASR; % the thing to plot
    zMax = max(M(:,:,:,clustNum),[],'all');
    zMin = min([ 0, min(M(:,:,:,clustNum),[],'all')]);

    % fig = myfig(0.4,'fig');
    fig = figure;%(3);

    nRows = 1;
    nNSets = 1;

    Colour = 'k';%jet(NFreq);
    % Colour = parula(NFreq);

    % xMin = -min([Stm(sel).PreT]) * 1e-3;
    % xMax = max([Stm(sel).StimT]+[Stm(sel).PostT]) * 1e-3;

    for s = 1:NSets
        setNum = 1; %plotSet(s);
        setIdx = 1; %find(FRA.FRASetNum == setNum);

        %FRA
        h = subplot(nRows, nNSets, s+0*NSets, 'Parent', fig);
        cla(h);

        % color map of spike rate
        CData = M(:, :, setIdx, clustNum);
        imagesc(h, CData, 'AlphaData', ~isnan(CData), [zMin,zMax]);

        % % contour of FACA p-value
        % Cont = -log10(FRA.FACApval(:, :, setIdx, clustNum));
        % hold(h,'on');
        % contour(h, Cont, [2, 3], 'w', 'ShowText', 'on');
        % hold(h, 'off');
        %
        % % contour of max neighbour correlation
        % Cont = (FRA.MaxNeighCorr(:, :, setIdx, clustNum));
        % hold(h,'on');
        % contour(h, Cont, [0.1, 0.2, 0.5], 'r', 'ShowText', 'on');
        % hold(h,'off');

        % format and label graph
        % title(h, [num2str(reTime(setIdx), '%.2f') ' h']);
        set(h,'Xscale', 'lin', 'YDir', 'normal',...
            'FontName', 'Arial', 'FontWeight', 'bold','FontSize', 12, ...
            'XTick',2:4:NFreq,'XTickLabel',round(UFreq(2:4:NFreq),1), 'XTickLabelRotation',45,...
            'YTick',2:2:NInt,'YTicklabel',UInt(2:2:NInt));
        xlabel(h,'Stimulus frequency (kHz)')

        if s==1
            ylabel(h, 'Stimulus intensity (dB SPL)')
        end
        cb = colorbar(h, 'eastoutside');
        cb.Label.String = 'Spike rate (Hz)';


    end
    sgtitle(['Cluster ' num2str(cids(clustNum))]) % whole figure title
    %if (length(spiketimes{clustNum}) < 500); close(gcf); end
end


%% PLOTFRASEQ plots a sequence of FRA
% Input:
%   FRA         = struct; containing results of analysis
%   clustNum    = scalar; cluster to plot (indexed 1:NClu)
%   RefTime     = [optional] reference time for display (e.g. Salicylate
%                 injection)
%   SpkT        = cell array of spike times relative to each. Complementary
%                 to Stm.
%   Stm         = struct array of stimulus parameters. To be made
%                 compatible with table format in the future.
%   plotSet     = vector containing the index of stimulus sets to be
%                 ploted. Indexed according the "SetNum".

axFRA = [];
axRast = [];
FRA = [];
clustNum = cids(:);
RefTime = [];
spiketimes = aligned_spikes;
Stm = stimuli_parameters.Stm;
plotSet = [];

plotFRASummary(axFRA,axRast,FRA, clustNum,RefTime,spiketimes,Stm, plotSet)

%% FRA analysis OLD
selSets = [];

% FRA_Stm = makeStmSets(S_sets,'FRA',selSets,ResPath);
%
% SpkT    =   AlignSpkS(FRA_Stm,SpkS,'FRA');
xf
FRA_Stm = 1;
spiketimes = aligned_spikes;

disp('Analyzing FRAs...'); tic
reOffset = [0,1];
[FRA,S_sets] = FRA_Analysis(FRA_Stm,spiketimes); %,0,0.015,10e-3,selSets,reOffset);
toc; disp('done.')
%%
% save FRA analysis
if exist([UserPath 'FRA_analysis\' DirName],'dir') == 0
    mkdir([UserPath 'FRA_analysis\' DirName]);
    mkdir([UserPath 'FRA_analysis\' DirName '\Figures']);
    mkdir([UserPath 'FRA_analysis\' DirName '\Figures\FRARaster\']);
    mkdir([UserPath 'FRA_analysis\' DirName '\Figures\CF-BF-Thr\']);
end
typelabel   =   'both';
save([UserPath 'FRA_analysis\' DirName '\' DirName '_FRA_' typelabel '_data.mat'],'cids','cpos','FRA','FRA_Stm','spiketimes');
% save([ResPath '\Results\' Mouse '_FRA_' typelabel '_data.mat'],'cids','cpos','FRA','FRA_Stm','SpkT');
disp('saved');


%% Plot FRA results OLD

% cluster = 2;
% plotSet = [6,8];
plotSet = FRA_Stm.setNum(strcmp(FRA_Stm.recType,'FRA'));
InjTime = FRA_Stm.startTime(S_sets.setNum==plotSet(1));
% InjIdx = checkInj(S_sets,plotSet);
% InjTime = FRA_Stm.startTime(InjIdx);
for cluster = 1:size(spiketimes,2)
    close all
    Name = [DirName ' -- Unit ' num2str(cids(cluster)) ' on Ch' num2str(cpos(cluster))];
    F = plotFRASeq(FRA, cluster,[],spiketimes,table2struct(FRA_Stm.Stm), plotSet);
    sgtitle(F,Name,'FontName','Arial','FontWeight','bold','FontSize', 16);

    if length(plotSet)==2
        set(F,'Position',[0,0,1600,1000]);
    elseif length(plotSet)==1
        set(F,'Position',[0,0,800,1000]);
    end

    saveas(F,[UserPath 'FRA_analysis\' DirName '\Figures\FRARaster\' DirName '-E' num2str(cpos(cluster)) '-U' num2str(cids(cluster)) '-FRA.png']);
end

%%local functions
function [sponR] = spontaneousrate(SpkS,S_sets,SetSamp)

USets = S_sets.setNum;
SetRec = S_sets.recType;
nSets = S_sets.NSets;
paths = S_sets.path;

nClust = size(SpkS,2);

sponR = nan(1,nClust);

File = paths{1,1};
Rec = string(File(1,end-6:end-4));
% Srise = Stm.Srise;
if Rec == "SIL"
    Samp = SetSamp(1);
elseif Rec == "FRA" || Rec == "AMn"
    load(paths{1,1},'Stm');
    Srise = Stm.Srise;
    Srise = Srise(1);
elseif Rec == "DRC" || Rec == "OMI"
    load(paths{1,1},'DigSamp');
    Srise = DigSamp{1,2}; Srise=Srise(1,1);
end

if Rec == "SIL"
    for s=1:nClust
        sponSpkS = SpkS{1,s};
        sponR(1,s) = (length(sponSpkS)/(Samp/30000));
    end
else
    for s=1:nClust
        %         S = SpkS{USets(k),s};
        S = SpkS{1,s}; %to make it work for chronic
        sponSpkS = S(S<=Srise); %extract all spikes prior to the start of the first trial of the set
        sponR(1,s) = (length(sponSpkS)/(Srise/30000)); %divide the spike count from this period by measurement period and store in sepSFR

    end
end


end

function [baseFR] = baselinerate(SpkS,Stm,PreT)

Set = Stm.Set;
USets = unique(Set);
nSets=length(USets);

nClust = size(SpkS,2);


baseFR = nan(nSets,nClust);

Srise = Stm.Srise;

for k=1:nSets
    sel = Set == USets(k);
    SetSrise = Srise(sel);
    nTrls = length(SetSrise);

    for s=1:nClust
        S = SpkS{USets(k),s};

        SCnt = zeros(nTrls,1);
        %calculating the spontaneous rate from the PreT period of trials
        for g=1:nTrls
            SCnt(g,1) = length(find(S <= SetSrise(g) & S >= SetSrise(g)-(PreT*30000)));
        end

        baseFR(k,s) = sum(SCnt)/(abs(PreT)*nTrls);



    end

end

end

function scnt   =   countSpks(SpkT,winStart,winEnd)
if (~iscell(SpkT)); SpkT = {SpkT}; end
scnt = nan(length(SpkT),1);
for t = 1:length(SpkT)
    if (isnan(SpkT{t}))
        continue
    end
    if(isempty(SpkT{t}))
        scnt(t) = 0;
    else
        scnt(t) = sum(SpkT{t} > winStart & SpkT{t} < winEnd);
    end
end
end

function Q = plotFRam(cluster, plotSet, SponR,BaseR,diffFR,ratFR,FRepoch,AMRes,Name)
Sets        =   AMRes.SetNum;
NSets       =   length(plotSet);
UMf         =	AMRes.UMf;
NMf         =	AMRes.NMf;
UMd         =	AMRes.UMd;
NMd         =	AMRes.NMd;
UInt        =   AMRes.UInt;
NInt        =   AMRes.NInt;

Col         =   [1 0 0; 0 0.5 0; 0 0 1; 0.75 0 0.75];


for s = 1:NSets
    myfig(1,'fig');
    cnt = 1;
    setSponR = BaseR(Sets==(plotSet(s)),cluster);
    for i = NInt:-1:1
        for d = 2:NMd
            mx      =   nan(NMf,1);
            subplot(NInt,NMd-1,cnt);
            plot([0 0],[-200 200],'-','Color',[0.4 0.4 0.4]);
            hold on;
            plot([-200 200],[0 0],'-','Color',[0.4 0.4 0.4]);
            X = FRepoch(UMd==0,UMf==0,UInt==(UInt(i)),Sets==(plotSet(s)),cluster,1)-setSponR;
            Y = FRepoch(UMd==0,UMf==0,UInt==(UInt(i)),Sets==(plotSet(s)),cluster,2)-setSponR;
            plot(X,Y,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8);
            mx(1) = max([abs(X) abs(Y)]);

            for f= 2:NMf
                X = FRepoch(UMd==(UMd(d)),UMf==(UMf(f)),UInt==(UInt(i)),Sets==(plotSet(s)),cluster,1)-setSponR;
                Y = FRepoch(UMd==(UMd(d)),UMf==(UMf(f)),UInt==(UInt(i)),Sets==(plotSet(s)),cluster,2)-setSponR;

                plot(X,Y,'o','MarkerEdgeColor',Col(f-1,:),'MarkerFaceColor',Col(f-1,:),'MarkerSize',8);
                mx(f) = max([abs(X) abs(Y)]);


            end

            mx  =   max(mx);
            ax = gca;
            plot([-200 200],[-200 200],'k--');
            set(ax,'YLim',[mx*-1.1 mx*1.1],'XLim',[mx*-1.1 mx*1.1]);
            xlabel('FRbase - FRum');
            ylabel('FRbase - FRam');
            title(['Int=' num2str(UInt(i)) ' Md=' num2str(UMd(d))]);

            cnt = cnt + 1;
        end
    end
    sgtitle([Name ' -- Set ' num2str(plotSet(s))]);
end

end

function minThr = FRAThreshold(setIdx,clusterIdx,FRA,FRA_Stm,Mouse,cpos)

deltaInt = [0,10,20,30,40];
%     Color = linspace(1,0,length(deltaInt))' * [0.5,0.5,1] ;%;0.25,0.25,.75;0,0,0];
Color =[0,      0,      1; ...
    0.25,   0.25,   1;...
    0.5,    0.5,    1; ...
    0.5,    0.5,    0.5;...
    0,      0,      0] ;%;0.25,0.25,.75;0,0,0];


for Set = setIdx
    setNum = FRA.FRAIdx(Set);
    for cluster = clusterIdx
        % set alpha value
        alpha = 1e-3;

        % assign significance (based on bootstrapped p-value)
        sig = double(FRA.FACApval(:,:,Set,cluster) <= alpha);
        sig(isnan(FRA.FACApval(:,:,Set,cluster))) = nan;

        % find minimum threshold and CF
        [row,col] = find(sig==1);
        if isempty(row)
            minThr = NaN;
        else
            minThr = FRA.UInt(min(row));
        end
    end
end
end

function F =  plotAMnRasterSel(AMRes,Name,T,AM_Stm,plotSet,RefTime)
NSets       =	length(plotSet);
NClu       	=	AMRes.NClu;

Stm         =   AM_Stm.Stm;
Mf          =   Stm.Mf;
Md          =   Stm.Md;
Int         =   Stm.Intensity;
Set         =   Stm.Set;



UMf         =	AMRes.UMf;
NMf         =	AMRes.NMf;
UMd         =	AMRes.UMd;
NMd         =	AMRes.NMd;
UInt        =   AMRes.UInt;
NInt        =   AMRes.NInt;

AMOnset     =   Stm.TransTime; AMOnset = unique(AMOnset);

startTime = AM_Stm.startTime;
meanTime = AM_Stm.meanTime;

if (isempty(RefTime)); RefTime = min(startTime);end
reTime = hours(meanTime - RefTime);


Colour      =   [1 0 0; 0 0.5 0; 0 0 1; 0.75 0 0.75];




F=myfig(0.8,'fig');

for b=1:NSets
    setNum = plotSet(b);
    setIdx = find(AM_Stm.setNum == setNum);

    subplot(1,NSets,b);
    cnt = 1;
    tick = 1;
    YTicks  =   nan((1+(1*(1)*(NMd-1))),1);
    YTickLab=   cell((1+(1*(1)*(NMd-1))),1);

    plotStore   =   cell((1+(1*(1)*(NMd-1))),1);
    MfLabIdx    =   zeros((1+(1*(1)*(NMd-1))),1);
    for k=1%:NInt
        sel = Md == 0 & Int == UInt(k) & Set == plotSet(b);
        nTrls = sum(sel);
        t  =   T(sel,1);

        for m=1:nTrls
            tt = t{m,1};
            if k==1 && m == 1
                MfLabIdx(tick,1) = 1;
            end

            plotStore{tick,1} = plot(tt,(cnt*ones(length(tt),1)),'k.','MarkerSize',4);
            hold on

            if m== nTrls/2
                YTicks(tick,1) =  cnt;
                YTickLab{tick,1} = ['UM'];
                tick = tick+1;
            end

            cnt = cnt+1;
        end

        cnt = cnt + 1;
    end
    cnt = cnt +2;

    for k=2
        Col = Colour(k-1,:);
        for m=1%:NInt
            for p=2:NMd
                sel = Md == UMd(p) & Mf == UMf(k) & Int == UInt(m) & Set == plotSet(b);
                nTrls = sum(sel);
                t   =   T(sel,1);

                for z=1:nTrls
                    if p==2 && z==1
                        MfLabIdx(tick,1) = 1;
                    end

                    tt = t{z,1};
                    plotStore{tick,1} = plot(tt,(cnt*ones(length(tt),1)),'.','MarkerEdgeColor', Col, 'MarkerFaceColor', Col,'MarkerSize',4);

                    if z== nTrls/2
                        YTicks(tick,1) =  cnt;
                        YTickLab{tick,1} = [num2str(UMd(p))];
                        tick = tick+1;
                    end

                    cnt = cnt+1;
                end

                %             keyboard
                cnt = cnt+2;
            end
            cnt = cnt+2;
        end
        cnt = cnt+2;
    end

    plot([AMOnset AMOnset],[0 cnt+5],'k-');
    a=gca;
    set(a,'FontSize',16,'FontName','Arial');
    ylim([0 cnt+3]);
    yticks(YTicks);
    yticklabels(YTickLab);
    xlim([-0.1 2.4]);
    %     title([num2str(reTime(setIdx),'%.2f') ' h']);



end

% sgtitle(Name,'FontName','Arial','FontSize',24);



if NMf == 1

elseif NMf == 2

elseif NMf == 3

elseif NMf == 4

elseif NMf == 5
    an = annotation('textbox',[0.2 0.5 0.3 0.3],'String',{['\color[rgb]{' num2str(Colour(end,:)) '}' num2str(UMf(end)) ' Hz'];...
        ['\color[rgb]{' num2str(Colour(end-1,:)) '}' num2str(UMf(end-1)) ' Hz'];...
        ['\color[rgb]{' num2str(Colour(end-2,:)) '}' num2str(UMf(end-2)) ' Hz'];
        ['\color[rgb]{' num2str(Colour(end-3,:)) '}' num2str(UMf(end-3)) ' Hz'];
        '\color[rgb]{0 0 0}UM'},'FitBoxToText','on','FontSize',14,'Position',[0.93, 0.8, 0.045052082145897,0.134994804178319],'EdgeColor','none');

end

end
function Wav = ChannelWaveforms2(KSpath,Mouse,ST,samp)

% ChMap       =   1+ [26    24    22    19    21    30    28    17    29    31    14    20    18    23    25    27     1     3     4 6     8    10    12    16    15    13    11     9     7     5     2     0    63    61    58    56    54    52 50    46    49    51    53    55    57    59    60    62    36    38    40    45    43    32    34    47    35 33    48    42    44    41    39    37];
% load('W:\KS\configfiles\CambridgeH3chronic.mat','chanMap');
load('W:\KS\configfiles\CambridgeH3acute.mat','chanMap');
ChMap = chanMap;
% ids     =   readNPY([KSpath '\spike_clusters.npy']);
% ST      =   readNPY([KSpath '\spike_times.npy']);
% amp     =   readNPY([KSpath '\amplitudes.npy']);
% if isempty(nWav)
nWav    =   length(ST);
% end

% sel     =   ids == cid;
% ST      =   ST(sel);
ids     =   ones(size(ST));
% amp     =   amp(sel);
% samp    =   -30:31;
t       =   (samp/30000)*1000;

gwfparams = struct;

gwfparams.dataDir = KSpath;    % KiloSort/Phy output folder
gwfparams.fileName = [Mouse '.dat'];         % .dat file containing the raw
%  gwfparams.dataDir = 'W:\temp';    % KiloSort/Phy output folder
%  gwfparams.fileName = 'temp_wh.dat';         % .dat file containing the raw
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 64;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [samp(1) samp(end)];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = nWav;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    ST; % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = ids; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

wf = getWaveForms(gwfparams,ChMap);


nSamp = length(samp);
mWav = wf.waveFormsMean; mWav = reshape(mWav,[64 nSamp]); mWav = 0.195*mWav;
Wav = wf.waveForms; Wav = reshape(Wav,[nWav,64,nSamp]); Wav = 0.195*Wav;
end
function F = ChannelWaveforms(F,KSpath,Mouse,cid,cpos,nWav)
% a function to load and plot raw waveforms.

% --- default parameter(s) ---
if isempty(nWav);    nWav    =   200; end
% ----------------------------

% --- Channel map ---
% ChMap       =   1+ [26    24    22    19    21    30    28    17    29    31    14    20    18    23    25    27     1     3     4 6     8    10    12    16    15    13    11     9     7     5     2     0    63    61    58    56    54    52 50    46    49    51    53    55    57    59    60    62    36    38    40    45    43    32    34    47    35 33    48    42    44    41    39    37];
% load('W:\KS\configfiles\CambridgeH3chronic.mat','chanMap');
load('W:\KS\configfiles\CambridgeH3acute.mat','chanMap');
ChMap = chanMap;
% -------------------

% --- Read info of individual (sorted) spikes ---
% all spikes
ids     =   readNPY([KSpath '\spike_clusters.npy']);
ST      =   readNPY([KSpath '\spike_times.npy']);
amp     =   readNPY([KSpath '\amplitudes.npy']);

% select spike with specified cluster ID "cid"
sel     =   ids == cid;
ST      =   ST(sel);
ids     =   ids(sel);
amp     =   amp(sel);
% -----------------------------------------------

% --- Parameters ---
samp    =   -30:31; % number of samples to plot (time);
NAdjCh  =   2; % number of adjacent channels (each side) to plot
sampRate = 30000;
t       =   (samp/sampRate)*1000;
nSamp   = length(samp);

gwfparams = struct;

gwfparams.dataDir = KSpath;    % KiloSort/Phy output folder
gwfparams.fileName = [Mouse '.dat'];         % .dat file containing the raw
%  gwfparams.dataDir = 'W:\temp';    % KiloSort/Phy output folder
%  gwfparams.fileName = 'temp_wh.dat';         % .dat file containing the raw
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 64;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [samp(1) samp(end)];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = nWav;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    ST; % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = ids; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

wf = getWaveForms(gwfparams,ChMap);

mWav = wf.waveFormsMean; mWav = reshape(mWav,[64 62]); mWav = 0.195*mWav;
Wav = wf.waveForms; Wav = reshape(Wav,[nWav,64,62]); Wav = 0.195*Wav;

%make subplot handles
%  w1 = subplot(4,4,1);
%  w2 = subplot(4,4,5);
%  w3 = subplot(4,4,9);
w1 = subplot(4,4,[1,9]);

a = subplot(4,4,13:16);
plot(a,(ST/30000),amp,'k.'); ylabel(a,'amplitude (a.u)'); xlabel(a,'time re: experiment');
set(a,'FontName','Arial','FontSize',12);
ylim(a,[0,Inf]);

W = subplot(4,4,[2.5:4 6.5:8 10.5:12]);

ch2plot = cpos + (-NAdjCh:NAdjCh);
ch2plot = ch2plot(ch2plot>0 & ch2plot <= 64);
mn = min(mWav(ch2plot,:),[],'all');
mx = max(mWav(ch2plot,:),[],'all');
offset = 1.2*(mx-mn);
plot(w1,t,mWav(ch2plot,:)' + ch2plot*offset,'k-')

%  if cpos == 1
%      plot(w1,t,mWav(cpos+1,:),'k-'); title(w1,['E' num2str(cpos+1)]);
%      plot(w2,t,mWav(cpos,:),'k-'); title(w2,['E' num2str(cpos)]);
%  elseif cpos == 64
%      plot(w2,t,mWav(cpos,:),'k-'); title(w2,['E' num2str(cpos)]);
%      plot(w3,t,mWav(cpos-1,:),'k-');title(w3,['E' num2str(cpos-1)]);
%  else
%      plot(w1,t,mWav(cpos+1,:),'k-'); title(w1,['E' num2str(cpos+1)]);
%      plot(w2,t,mWav(cpos,:),'k-'); title(w2,['E' num2str(cpos)]);
%      plot(w3,t,mWav(cpos-1,:),'k-');title(w3,['E' num2str(cpos-1)]);
%  end

%  ylabel(w2,'Vpp (\muV)');
% yticks(w1,ch2plot*offset);yticklabels(ch2plot);
set(w1,'YTick',ch2plot*offset,'YTickLabel',ch2plot)
xlabel(w1,'time (ms)');
set(w1,'FontName','Arial','FontSize',12);
%  set(w2,'FontName','Arial','FontSize',12);
%  set(w3,'FontName','Arial','FontSize',12);
%  linkaxes([w1 w2 w3],'y');
%  mn = min(mWav(cpos,:));
%  mx = max(mWav(cpos,:));
%  df = diff([mn mx]);
%  ylim(w2,[mn-(df/10) mx+(df/10)]);

Wo = reshape(Wav(:,cpos,:),[nWav, nSamp]);
for k=1:nWav
    w=Wo(k,:) - (mean(Wo(k,1:20)));
    plot(t,w,'k-','LineWidth',0.5);
    hold on;
end
xlim(W,[min(t) max(t)]);
set(W,'FontName','Arial','FontSize',12);
ylabel('Vpp (\muV)');
xlabel('time (ms)');

sgtitle([Mouse ' - unit ' num2str(cid)],'FontName','Arial','FontSize',25);
end
%% Local functions % OLD


function plotFRASummary(axFRA,axRast,FRA, clustNum,RefTime,SpkT,Stm, plotSet)
%% PLOTFRASEQ plots a sequence of FRA
% Input:
%   FRA         = struct; containing results of analysis
%   clustNum    = scalar; cluster to plot (indexed 1:NClu)
%   RefTime     = [optional] reference time for display (e.g. Salicylate
%                 injection)
%   SpkT        = cell array of spike times relative to each. Complementary
%                 to Stm.
%   Stm         = struct array of stimulus parameters. To be made
%                 compatible with table format in the future.
%   plotSet     = vector containing the index of stimulus sets to be
%                 ploted. Indexed according the "SetNum".

if nargin< 6 || isempty(plotSet); plotSet = 1; end %FRA.FRASetNum;end

NSets       =	length(plotSet);
NClu       =	FRA.NClu;
UFreq       =	FRA.UFreq;
NFreq       =	FRA.NFreq;
UInt       =	FRA.UInt;
NInt       =	FRA.NInt;

startTime = FRA.startTime;
meanTime = FRA.meanTime;

M   =   FRA.FRASR; % the thing to plot
zMax = max(M(:,:,:,clustNum),[],'all');
zMin = min([ 0, min(M(:,:,:,clustNum),[],'all')]);

if (nargin < 3 || isempty(RefTime)); RefTime = min(startTime);end
reTime = hours(meanTime - RefTime);

% fig = myfig(0.4,'fig');

% fig = figure;%(3);



nRows = 2;

Colour = 'k';%jet(NFreq);
% Colour = parula(NFreq);

% check duration
sel = zeros(size(Stm));
%
%     tempSel = reshape([Stm.Set] == plotSet(s),size(sel));
%     sel = sel | tempSel;
sel = Stm.Set == plotSet;

xMin = -min([Stm.PreT(sel)]) *(1e-3);
xMax = max([Stm.StimT(sel)]+[Stm.PostT(sel)]) * 1e-3;

% xMin = -min([Stm(sel).PreT]) * 1e-3;
% xMax = max([Stm(sel).StimT]+[Stm(sel).PostT]) * 1e-3;

for s = 1:NSets
    setNum = plotSet(s);
    setIdx = find(FRA.FRASetNum == setNum);

    %FRA
    h = axFRA;
    % color map of spike rate
    CData = M(:,:,setIdx,clustNum);
    imagesc(h,CData,'AlphaData',~isnan(CData),[zMin,zMax]);
    % contour of FACA p-value
    Cont = -log10(FRA.FACApval(:,:,setIdx,clustNum));
    hold(h,'on');contour(h,Cont,[2,3],'w','ShowText','on');hold(h,'off');
    % contour of max neighbour correlation
    Cont = (FRA.MaxNeighCorr(:,:,setIdx,clustNum));
    hold(h,'on');contour(h,Cont,[0.1,0.2,0.5],'r','ShowText','on');hold(h,'off');
    % format and label graph
    %     title(h,[num2str(reTime(setIdx),'%.2f') ' h']);
    set(h,'Xscale','lin','YDir','normal',...
        'XTick',2:4:NFreq,'XTickLabel',round(UFreq(2:4:NFreq),1), 'XTickLabelRotation',45,...
        'YTick',2:2:NInt,'YTicklabel',UInt(2:2:NInt));
    xlabel(h,'frequency (kHz)')
    if s==1
        ylabel(h,'intensity (dB SPL)')
    end
    cb = colorbar(h,'eastoutside');
    cb.Label.String = 'Spike Rate (Hz)';

    ax = gca;
    set(ax,'FontSize',10,'FontName','Arial');
    %Raster Plot

    h = axRast;
    cla(h);
    sel = [Stm.Set] == setNum;
    [h,YTick,YTickLab] = plotraster(h, SpkT(sel,clustNum),[[Stm.Freq(sel)], [Stm.Intensity(sel)]], Colour, [100;0],1);
    yticks = YTick{1};
    set(h,'YTick',yticks(2:4:end),'YTicklabel',round(UFreq(2:4:end),1),...
        'XTick',0:0.05:xMax);
    ax = gca;
    set(ax,'FontSize',10,'FontName','Arial');
    %    set(h,'YTick',yticks([1,2:2:end]),'YTicklabel',round(UFreq([1,2:2:end]),1));
    xlim(h,[xMin,xMax]);
    margin = 100;
    ylim(h,[min(YTick{end})-margin max(YTick{end})+margin]);
    xlabel('time re: stimulus (sec)');
    if s==1
        ylabel(h,'frequency (kHz)');
    end
end

% disp(zMin);
% disp(zMax);
end