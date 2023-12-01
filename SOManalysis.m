function SOM = SOManalysis(stimuli_parameters, aligned_spikes, cids)
%% SOM analysis
% input: stimuli_parameters.Par, stimuli_parameters.Stm, aligned_spikes
% output: first spike latency SOM/AM trials

% IN PROGRESS

% stimulus parameters
NClu = length(cids); % cluster info
UAmp = unique([stimuli_parameters.Stm.Amplitude]);
NAmp = length(UAmp);

% analysis period
winStart = 0;
winEnd   = 0.115;

% % initiation variables
%FRAScntSD = nan(NInt, NFreq, NSets, NClu); % spike count
%FRAScnt = nan(NInt, NFreq, NSets, NClu); % spike count
%FRASR = nan(NInt, NFreq, NSets, NClu); % spike count
NTrials = nan(NAmp, NClu); % number of trials
MedFSL = nan(NAmp, NClu); % median first spike latency

for cluster = 1:NClu
    disp(['Analysing cluster ' num2str(cluster) ' of ' num2str(NClu)]);

    % Spike count analysis
    for amplitude = 1:NAmp
            sel = [stimuli_parameters.Stm.Amplitude] == UAmp(amplitude); % selects correct trials from params
            NTrials(amplitude,cluster) = sum(sel);

            if (NTrials(amplitude,cluster) == 0); continue; end

            %count spikes
            SCnt = nan(NTrials(amplitude,cluster), 1);

            tempSpiketimes = aligned_spikes(sel,cluster); % adjust for previous sessions, selects incorrect rows from alined_spikes
            for t = 1:length(tempSpiketimes)
                if (isnan(tempSpiketimes{t})); continue; end
                if(isempty(tempSpiketimes{t}))
                    SCnt(t) = 0;
                else
                    SCnt(t) = sum(tempSpiketimes{t} > winStart & tempSpiketimes{t} < winEnd);
                end
            end

            %FRAScntSD(intensity,freq,setNum,cluster)  = std(SCnt,'omitnan');
            %FRAScnt(intensity,freq,setNum,cluster)    = mean(SCnt,'omitnan');
            %FRASR(intensity,freq,setNum,cluster)      = mean(SCnt,'omitnan') / (winEnd - winStart);

            %FSL
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
            
            MedFSL(amplitude,cluster) = median(fsl);

    end % amplitude loop
    clearvars('tempSpiketimes');

end % cluster loop


% Output data
% Experiment meta data
SOM.NSets        = NSets;

% cluster
SOM.NClu = NClu;

% stimulus parameters
SOM.UFreq = UAmp;
SOM.NFreq = NAmp;
SOM.UInt = UInt;
SOM.NInt = NInt;

% results - number of trials
SOM.NTrials = NTrials;

% results - spike latency
SOM.MedFSL = MedFSL;

% results - spike count
SOM.FRASR = FRASR;
SOM.FRAScnt = FRAScnt;
SOM.FRAScntSD = FRAScntSD;

%% MedFSL heatmap
for clustNum = 1:NClu
    NSets       =	1;
    NClu       =	SOM.NClu;
    UAmp       =	SOM.UFreq;
    NAmp       =	SOM.NFreq;
    UInt       =	SOM.UInt;
    NInt       =	SOM.NInt;

    startTime = 0;%FRA.startTime;
    meanTime = 0;%FRA.meanTime;

    M   =   MedFSL; % the thing to plot
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

        % format and label graph
        % title(h, [num2str(reTime(setIdx), '%.2f') ' h']);
        set(h,'Xscale', 'lin', 'YDir', 'normal',...
            'FontName', 'Arial', 'FontWeight', 'bold','FontSize', 12, ...
            'XTick',2:4:NAmp,'XTickLabel',round(UAmp(2:4:NAmp),1), 'XTickLabelRotation',45,...
            'YTick',2:2:NInt,'YTicklabel',UInt(2:2:NInt));
        xlabel(h,'Stimulus frequency (kHz)')

        if s==1
            ylabel(h, 'Stimulus intensity (dB SPL)')
        end
        cb = colorbar(h, 'eastoutside');
        cb.Label.String = 'Spike rate (Hz)';

    end
    sgtitle(['FSL: Cluster ' num2str(cids(clustNum))]) % whole figure title
    %if (length(spiketimes{clustNum}) < 500); close(gcf); end
end
