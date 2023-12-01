function FRA = FRAanalysis(stimuli_parameters, aligned_spikes, cids, OutPath)
% FRA analysis
% input: stimuli_parameters.Par, stimuli_parameters.Stm, aligned_spikes
% output: FRA and first spike latency

NClu = length(cids); % cluster info
NSets = 1;
setNum = 1; % set info

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

            tempSpiketimes = aligned_spikes(sel,cluster); % spike times per combination per cluster
            for t = 1:length(tempSpiketimes)
                if (isnan(tempSpiketimes{t})); continue; end
                if(isempty(tempSpiketimes{t}))
                    SCnt(t) = 0;
                else
                    SCnt(t) = sum(tempSpiketimes{t} > winStart & tempSpiketimes{t} < winEnd);
                end
            end

            FRAScntSD(intensity,freq,setNum,cluster)  = std(SCnt,'omitnan');
            FRAScnt(intensity,freq,setNum,cluster)    = mean(SCnt,'omitnan');
            FRASR(intensity,freq,setNum,cluster)      = mean(SCnt,'omitnan') / (winEnd - winStart);

            %FSL
            fsl = inf(length(tempSpiketimes), 1); % intiate
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
            
            MedFSL(intensity, freq, setNum, cluster) = median(fsl);

        end % intensity loop
    end % frequency loop
    clearvars('tempSpiketimes');

end % cluster loop

% Output data
% Experiment meta data
FRA.NSets        = NSets;

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
    if zMax == 0; continue; end
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

    sgtitle(['FRA (unit ' num2str(cids(clustNum)) ')']) % whole figure title
    %if (length(spiketimes{clustNum}) < 500); close(gcf); end

    figname = sprintf('M02_FRA_cluster %i', cids(clustNum));
    saveas(gcf, fullfile(OutPath, [figname '.jpg']));
    saveas(fig, fullfile(OutPath, figname));

end


%% MedFSL heatmap
for clustNum = 1:NClu
    NSets       =	1;
    NClu       =	FRA.NClu;
    UFreq       =	FRA.UFreq;
    NFreq       =	FRA.NFreq;
    UInt       =	FRA.UInt;
    NInt       =	FRA.NInt;

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
        cb.Label.String = 'Latency first spike (s)';



    end

    sgtitle(['FSL (unit ' num2str(cids(clustNum)) ')']) % whole figure title
    figname = sprintf('M02_FRA FSL_cluster %i', cids(clustNum));
    saveas(gcf, fullfile(OutPath, [figname '.jpg']));
    saveas(fig, fullfile(OutPath, figname));

    %if (length(spiketimes{clustNum}) < 500); close(gcf); end
end
